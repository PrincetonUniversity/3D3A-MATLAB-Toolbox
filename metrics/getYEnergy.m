function [DE,CDE,BX] = getYEnergy(YC,varargin)
%GETYENERGY Compute the energy of spherical harmonic coefficients.
%   [DE,CDE,BX] = GETYENERGY(YC) takes a matrix of spherical harmonic
%   coefficients, YC, and returns directional energy, DE, and cumulative
%   directional energy, CDE, spectra at and up to, respectively, each 
%   spherical harmonic degree. DE and CDE are specified in dB.
%
%   If YC is a matrix with dimensions (DEG+1)^2-by-N, with N >= 1, then DE
%   and CDE will each have dimensions (DEG+1)-by-N. Also returned is BX, a 
%   length-N vector of the smallest spherical harmonic degree at each 
%   frequency where the value of CDE is greater than or equal to X = 99% 
%   of the total energy across all degrees.
%
%   [DE,CDE,BX] = GETYENERGY(YC,X) optionally allows the specification of X
%   used to compute BX (see above). X must be a real scalar between 0 and
%   100. If X is specified as [], a default value of 99 is used.
%
%   [DE,CDE,BX] = GETYENERGY(YC,X,NORMVAL) optionally specifies a 
%   normalization value to which the maximum of all coefficients (at each 
%   frequency) is normalized. Not specifying a value results in un-
%   normalized values. NORMVAL should be specified in dB.
%
%   [DE,CDE,BX] = GETYENERGY(YC,X,NORMVAL,NORMTYPE) optionally specifies
%   the type of normalization to perform. The two options are:
%       1. 'max' (default) - The value of DE at each frequency is
%       normalized by the maximum value of DE at that frequency. This type
%       of normalization allows us to easily see, for a given frequency, 
%       which degree contains the most energy. The relative total energy
%       across frequencies remains unchanged (effectively, the energy
%       spectrum for any given source direction is not "flattened").
%
%       2. 'sum' - The value of DE at each frequency is normalized by the
%       sum of all DE at that frequency. This type of normalization allows 
%       us to easily see, for a given frequency, the distribution of energy
%       across degree, irrespective of the total energy at that frequency.
%       The relative total energy across frequencies is modified (
%       effectively, the energy spectrum for any given source direction is 
%       "flattened").
%
%   In both cases, CDE is computed directly from the normalized DE values.
%   This option is relevant only when a value for NORMVAL is specified.

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

%   References:
%       [1] Brinkmann and Weinzierl (2018) - Comparison of head-related 
%   transfer functions pre-processing techniques for spherical harmonics 
%   decomposition.

%% Check inputs

narginchk(1,4)

if nargin < 4
    NORMTYPE = 'max';
else
    NORMTYPE = varargin{3};
    validateattributes(NORMTYPE,{'char'},{'scalartext','nonempty'},...
        'getYEnergy','NORMTYPE',4)
end

if nargin < 3
    NORMVAL = []; % Indicates not to perform normalization
else
    NORMVAL = varargin{2};
    validateattributes(NORMVAL,{'numeric'},{'scalar','nonnan','finite',...
        'real'},'getYEnergy','NORMVAL',3)
end

if nargin < 2
    X = 99;
else
    X = varargin{1};
    if isempty(X)
        X = 99;
    else
        validateattributes(X,{'numeric'},{'scalar','nonnan','finite',...
        'real','nonnegative','<=',100},'getYEnergy','X',2)
    end
end

[nR,nC] = size(YC);
DEG = round(sqrt(nR)-1);
if (DEG+1)^2 ~= nR
    error(['Invalid YC data. The number of rows in YC must be',...
        ' (DEG+1)^2 for a non-negative integer, DEG.'])
end

%% Compute DE and CDE

degVec = 0:DEG;
YCEnergies = abs(YC).^2;
DE = zeros(DEG+1,nC);
for ii = 0:DEG
    lN = (degVec(ii+1)^2)+1;
    uN = (degVec(ii+1)+1)^2;
    DE(ii+1,:) = sum(YCEnergies(lN:uN,:),1);
end

if ~isempty(NORMVAL)
    switch lower(NORMTYPE)
        case 'max'
            maxVals = repmat(max(DE,[],1),DEG+1,1);
        case 'sum'
            maxVals = repmat(sum(DE,1),DEG+1,1);
        otherwise
            error('Invalid NORMTYPE specification.')
    end
    
    DE = (DE./maxVals)*db2mag(NORMVAL);
end

CDE = cumsum(DE,1);

%% Compute BX

BX = zeros(nC,1);
for ii = 1:nC
    indx = find(CDE(:,ii) >= (X/100)*max(CDE(:,ii),[],1),1);
    BX(ii) = indx-1;
end

%% Convert to dB and return

DE = 10*log10(DE);
CDE = 10*log10(CDE);

end
