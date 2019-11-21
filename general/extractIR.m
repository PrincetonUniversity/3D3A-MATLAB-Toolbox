function [irOut,posOut,indx] = extractIR(irMat,posMat,posOut,varargin)
%EXTRACTIR Find nearest IRs to desired positions.
%   HOUT = EXTRACTIR(HIN,RIN,ROUT) given a set of IRs HIN with
%   corresponding positions RIN, this function returns HOUT, the subset of
%   IRs which are nearest to the desired positions ROUT. HIN should be a
%   matrix of size K-by-P, where K is the impulse response length and P is
%   the number of positions. Correspondingly, RIN should be a matrix of
%   size P-by-3. Then, if ROUT is size Q-by-3, HOUT will be K-by-Q. All
%   position vectors should be given in Cartesian coordinates.
%
%   HOUT = EXTRACTIR(HIN,RIN,ROUT,NORMFLAG) additionally specifies whether
%   the position vectors should be normalized to unit length. By default,
%   NORMFLAG is false.
%
%   HOUT = EXTRACTIR(HIN,RIN,ROUT,NORMFLAG,EXACTFLAG) additionally 
%   specifies whether ROUT is a subset of RIN or not. If EXACTFLAG is true, 
%   then an error is thrown when any set of coordinates in ROUT is not 
%   found in RIN. Otherwise, the function performs as expected.
%
%   HOUT = EXTRACTIR(HIN,RIN,ROUT,NORMFLAG,true,N) additionally specifies a 
%   precision value, N, to use when rounding coordinates in RIN and ROUT 
%   when EXACTFLAG is set to true. If N is specified, then the values in 
%   RIN and ROUT are first rounded to N digits before comparisons are made. 
%   If N is not specified and EXACTFLAG is set to true, then N is assumed 
%   to be 3. To use full precision (i.e., no rounding), set N to inf.
%
%   [HOUT,ROUT,INDX] = EXTRACTIR(___) additionally returns the actual
%   nearest positions ROUT which exist in RIN as well as a column vector,
%   INDX, of the corresponding row-indices of RIN. By convention, INDX is
%   given such that HOUT = HIN(:,INDX) and ROUT = RIN(INDX,:).
%
%   See also FINDNEAREST, NORMALIZEVECTOR.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2019 Princeton University
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

% Check number of inputs
narginchk(3,6);

% Validate attributes of required inputs
validateattributes(irMat,{'numeric'},{'2d','nonempty'},'extractIR',...
    'HIN',1);
irMat = shiftdim(irMat); % If row vector, force to column.
numIRs = size(irMat,2);
validateattributes(posMat,{'numeric'},{'2d','nonempty','real','size',...
    [numIRs,3]},'extractIR','RIN',2);
validateattributes(posOut,{'numeric'},{'2d','nonempty','real','size',...
    [NaN,3]},'extractIR','ROUT',3);

% Parse optional inputs
switch nargin
    case 3
        presVal = 3;
        eFlag = false;
        normFlag = false;
    case 4
        presVal = 3;
        eFlag = false;
        normFlag = varargin{1};
    case 5
        presVal = 3;
        eFlag = varargin{2};
        normFlag = varargin{1};
    case 6
        presVal = varargin{3};
        eFlag = varargin{2};
        normFlag = varargin{1};
end
    
if isempty(eFlag)
    eFlag = false;
end

if isempty(normFlag)
    normFlag = false;
end

% Validate attributes of optional inputs
validateattributes(normFlag,{'logical'},{'scalar','nonempty'},...
    'extractIR','NORMFLAG',4);
validateattributes(eFlag,{'logical'},{'scalar','nonempty'},...
    'extractIR','EXACTFLAG',5);
validateattributes(presVal,{'numeric'},{'scalar','nonempty','nonnan'},...
    'extractIR','N',6);

if normFlag
    posMat = normalizeVector(posMat,2);
    posOut = normalizeVector(posOut,2);
end

numOuts = size(posOut,1);
indx = zeros(numOuts,1);
if eFlag
    if isfinite(presVal)
        posMat = round(posMat,presVal);
        posOut = round(posOut,presVal);
    end
    for ii = 1:numOuts
        indx(ii) = find(posMat(:,1) == posOut(ii,1) & ...
            posMat(:,2) == posOut(ii,2) & posMat(:,3) == posOut(ii,3),1);
        if isempty(indx(ii))
            error('Desired source position could not be found.')
        end
    end
else
    for ii = 1:numOuts
        [~,indx(ii)] = findNearest(posMat,posOut(ii,:),'l2');
    end
end

irOut = irMat(:,indx);
posOut = posMat(indx,:);

end
