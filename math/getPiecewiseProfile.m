function profileOut = getPiecewiseProfile(params,pLen)
%GETPIECEWISEPROFILE Compute a piecewise regularization profile.
%   profileOut = GETPIECEWISEPROFILE(params,pLen) returns a piecewise 
%   regularization profile of length pLen given input parameters, params, 
%   as an N-by-3 matrix, where N is the number of "pieces". Each row of 
%   params is specified as [w1,w2,eps], where w1 and w2 specify the 
%   approximate start and end points of a given piece and eps specifies the 
%   value/height of the piece. w1 and w2 must range from 0 to 1 with w2 > 
%   w1. The values of w1 and w2 in row i must exceed the corresponding 
%   values in row j, where i > j. eps must be non-negative. The profile is 
%   returned as a column vector of length pLen.

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
%   Copyright (c) 2018 Princeton University
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

narginchk(2,2);

% Validate inputs
validateattributes(params,{'double'},{'2d','nonempty','nonnan','finite',...
    'real','ncols',3},'getPiecewiseProfile','params',1)
validateattributes(pLen,{'double'},{'scalar','nonempty','nonnan',...
    'finite','integer','positive'},'getPiecewiseProfile','pLen',2)

evenFlag = false;
if floor(pLen/2) == (pLen/2) % If pLen is even
    evenFlag = true;
    pLen = pLen + 1;
end

% Check if w1 and w2 fall between 0 and 1
validateattributes(params(:,1),{'double'},{'increasing','nonnegative',...
    '<=',1},'getPiecewiseProfile','all values of w1 in params')
validateattributes(params(:,2),{'double'},{'increasing','nonnegative',...
    '<=',1},'getPiecewiseProfile','all values of w2 in params')

lIndxs = floor(params(:,1)*(pLen-1))+1;
hIndxs = floor(params(:,2)*(pLen-1))+1;
indxDiff = lIndxs(2:end)-hIndxs(1:(end-1));

% Check if indxDiff is non-decreasing
validateattributes(indxDiff,{'double'},{'nondecreasing'},...
    'getPiecewiseProfile','all values of w1-w2 in params')

profileOut = zeros(pLen,1); % Initialize output profile
numPieces = size(params,1);
for ii = 1:(numPieces-1)
    profileOut(lIndxs(ii):hIndxs(ii)) = params(ii,3);
    if lIndxs(ii+1) > hIndxs(ii)
        xRange = lIndxs(ii+1)-hIndxs(ii);
        xVec = (hIndxs(ii):lIndxs(ii+1)).';
        m = (params(ii+1,3)-params(ii,3))/xRange;
        profileOut(hIndxs(ii):lIndxs(ii+1)) = profileOut(hIndxs(ii))+m*...
            (xVec-hIndxs(ii));
    end    
end
profileOut(lIndxs(numPieces):hIndxs(numPieces)) = params(numPieces,3);

if evenFlag
    tsin = timeseries(profileOut,linspace(0,1,pLen));
    tsout = resample(tsin,linspace(0,1,pLen-1));
    profileOut = tsout.Data;
end

end
