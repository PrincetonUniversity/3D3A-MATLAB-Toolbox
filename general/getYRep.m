function [OD,YC,YM,DEG] = getYRep(ID,YD,varargin)
%GETYREP Compute a spherical harmonic representation of input data.
%   [OD,YC,YM,DEG] = GETYREP(ID,{'Ys',YM}) returns the spherical harmonic 
%   representation of the data in ID using the matrix of spherical 
%   harmonics in YM. If ID has dimensions N-by-p, YM must have dimensions 
%   p-by-M. See COMPUTEYMAT to compute YM. 
%
%   The input data reconstructed using spherical harmonics will be returned 
%   as OD, which will have the same dimensions as ID. Also returned are the 
%   spherical harmonic coefficients, YC, which will have dimensions M-by-N.
%   The coefficients in YC are computed using a least-squares approach
%   without any regularization. The spherical harmonic degree used to 
%   represent the input data is returned as DEG (especially useful when the 
%   'checkMaxD' option below is set to true). The input YM is also returned 
%   unaltered in this case.
%
%   ___ = GETYREP(ID,{'Ys',YM,'Wm',WV}) returns the spherical harmonic 
%   representation of the data in ID using the matrix of spherical 
%   harmonics in YM, and a vector of weights, WV. When WV is not specified, 
%   YC is computed using a least-squares minimization. Here, YC is computed 
%   using YM and WV directly (i.e., without using a least-squares 
%   approach). For example, if YM is specified on a quadrature grid with WV 
%   containing the quadrature weights, YC are obtained from these 
%   quantities directly rather than as the outcome of an optimization 
%   problem.
%
%   ___ = GETYREP(ID,{'YDirs',SP,'D',maxD}) returns the maxD-degree 
%   spherical harmonic representation of the data in ID whose values are 
%   specified for the spatial directions given in SP. If ID has dimensions 
%   N-by-p, SP must be specified in SOFA cartesian coordinates as a p-by-3 
%   matrix, where p is the number of directions. maxD must be a non-
%   negative scalar. In addition to returning OD and YC, the matrix of 
%   spherical harmonics, YM, computed internally using SP and maxD, are 
%   also returned. YM will have dimensions p-by-(maxD+1)^2. The following 
%   optional inputs may also be provided:
%       1. GETYREP(...,'TYPE','real') - use real-valued spherical
%       harmonics (default).
%       2. GETYREP(...,'TYPE','complex') - use complex-valued spherical
%       harmonics.
%       3. GETYREP(...,'CSPHASE',0) - Ignore Condon-Shortley phase
%       (default).
%       4. GETYREP(...,'CSPHASE',1) - Include Condon-Shortley phase.
%
%   ___ = GETYREP(___,'checkMaxD',VAL) - specify whether or not to limit 
%   the max. spherical harmonic degree to minimize the possibility of 
%   aliasing and truncation errors being introduced. VAL can take the
%   following values:
%       1. true/false - Do/do not limit degree. Use specified maxD or value
%       inferred from number of columns in YM. The default is false.
%
%       2. {true,L} - Limit degree using an oversampling factor, L >= 1.
%       Choosing this option may result in the specified value of maxD
%       being overridden. If L is not specified, a default value of 4 is
%       chosen. Some common values of L are:
%           i.      L = 1.3 - Data sampled on a Lebedev grid.
%           ii.     L = 2 - Data sampled using a Gaussian scheme.
%           iii.    L = 4 - Data sampled on an equiangular grid (i.e.,
%           uniform spacing in azimuth and elevation, and equal number of
%           azimuths and elevations).
%
%       For more information on oversampling factors, see Rafaely [1].
%
%   See also COMPUTEYMAT, COMPUTEY.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2021 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a
%   copy of this software and associated documentation files (the 
%   "Software"), to deal in the Software without restriction, including 
%   without limitation the rights to use, copy, modify, merge, publish, 
%   distribute, sublicense, and/or sell copies of the Software, and to 
%   permit persons to whom the Software is furnished to do so, subject to 
%   the following conditions:
%   
%   The above copyright notice and this permission notice shall be included
%   in all copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
%   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
%   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
%   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%   =======================================================================

%   Refs:
%       [1]. Boaz Rafaely (2005) - Analysis and Design of Spherical 
%   Microphone Arrays.

%% Check inputs

narginchk(2,8)

validateattributes(ID,{'numeric'},{'2d','nonempty','nonnan','finite'},...
    'getYRep','ID',1)
validateattributes(YD,{'cell'},{'nonempty'},'getYRep','YD',2)
numVarArgIn = length(varargin);

indx = find(strcmpi(varargin,'checkMaxD'),1);
if isempty(indx)
    checkMaxDFlag = false;
else
    if numVarArgIn >= indx+1
        VAL = varargin{indx+1};
        if iscell(VAL)
            checkMaxDFlag = VAL{1};
        else
            checkMaxDFlag = VAL;
        end
    else
        error('Option ''checkMaxDFlag'' must be followed by VAL.')
    end
    validateattributes(checkMaxDFlag,{'logical'},{'scalar'},'getYRep',...
        'checkMaxD')
end

optFlag = true;
[dataLen,numDirs] = size(ID);
indx = find(strcmpi(YD,'Ys'),1);
if ~isempty(indx)
    switch length(YD)
        case 2
            YM = YD{indx+1};
        case 4
            YM = YD{indx+1};
            indx2 = find(strcmpi(YD,'Wm'),1);
            if ~isempty(indx2)
                WV = shiftdim(YD{indx2+1});
            else
                error('Unrecognized entry in second input.')
            end
            optFlag = false;
        otherwise
            error('Invalid second input. Must contain 2 or 4 elements.')
    end
    maxD = sqrt(size(YM,2))-1;
else
    indx3 = find(strcmpi(YD,'YDirs'),1);
    if ~isempty(indx3)
        switch length(YD)
            case 4
                SP = YD{indx3+1};
                if strcmpi(YD{indx3+2},'D')
                    maxD = YD{indx3+3};
                else
                    error('Unrecognized entry in second input.')
                end
            otherwise
                error(['Invalid second input. Must contain 4 elements',...
                    ' when ''YDirs'' is specified.'])
        end
        
        indx4 = find(strcmpi(varargin,'TYPE'),1);
        if ~isempty(indx4)
            if numVarArgIn > indx4
                typeVal = varargin{indx4+1};
            else
                error('Option ''TYPE'' must be followed by a string.')
            end
        else
            typeVal = 'real';
        end
        
        indx5 = find(strcmpi(varargin,'CSPHASE'),1);
        if ~isempty(indx5)
            if numVarArgIn > indx5
                csPhaseFlag = varargin{indx4+1};
            else
                error('Option ''CSPHASE'' must be followed by a scalar.')
            end
        else
            csPhaseFlag = 0;
        end
    else
        error('Invalid second input. Must contain ''Ys'' or ''YDirs''.')
    end
end

%% Prepare for main calculation

if checkMaxDFlag
    if length(VAL) == 2
        L = VAL{2};
        validateattributes(L,{'numeric'},{'scalar','real','>=',1},...
            'getYRep','L')
    else
        L = 4;
    end
    DEG = min([maxD,floor(sqrt(numDirs/L))-1]);
else
    DEG = maxD;
end

if isempty(indx)
    % Compute matrix of spherical harmonics up to degree DEG sampled at 
    % directions specified in SP.
    YM = computeYMat(DEG,SP,'TYPE',typeVal,'CSPHASE',csPhaseFlag);
end

%% Perform main calculation

maxPickIndx = (DEG+1)^2;
pickIndxs = 1:maxPickIndx;
YM = YM(:,pickIndxs);
if optFlag
    invYMat = pinv(YM);
    YC = invYMat*(ID.');
else
    YC = 4*pi*(repmat(WV.',dataLen,1).*ID*conj(YM)).';
end
OD = (YM*YC).';

end
