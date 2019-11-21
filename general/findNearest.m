function [xOut,IND] = findNearest(xIn,xD,varargin)
%FINDNEAREST Find the row of a matrix nearest to a given row vector.
%   Y = FINDNEAREST(X,XD) returns Y, the row of the matrix X which is
%   closest (in a least-squares sense) to the row vector XD. The number of 
%   columns in X must equal the length of XD. If X is a vector, XD must be
%   a scalar.
%
%   Y = FINDNEAREST(X,XD,SCALE) optionally uses SCALE to measure the
%   distance between elements. The options for SCALE are:
%       1. 'linear' or 'l1' - distance is measured using the manhattan
%       norm.
%       2. 'least-squares' or 'l2' - distance is measured using the
%       Euclidean norm. This is the default.
%       3. 'angular' - distance is measured in terms of the angle between
%       two positions.
%
%   Y = FINDNEAREST(X,XD,SCALE,N) optionally returns the N nearest rows of
%   X.
%
%   [Y,I] = FINDNEAREST(___) also returns the position(s), I, of Y in X.

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
narginchk(2,4);

% Validate attributes of required inputs
validateattributes(xIn,{'numeric'},{'2d','nonempty','nonnan'},...
    'findNearest','X',1);
xIn = shiftdim(xIn); % If vector, force to a column.
[numRowsX,numColsX] = size(xIn);
validateattributes(xD,{'numeric'},{'row','nonempty','ncols',numColsX},...
    'findNearest','XD',2);

% Parse optional inputs
switch nargin
    case 2
        SCALE = 'l2';
        N = 1;
    case 3
        SCALE = varargin{1};
        N = 1;
    case 4
        SCALE = varargin{1};
        N = varargin{2};
end

if isempty(SCALE)
    SCALE = 'l2';
end

% Validate attributes of optional inputs
validateattributes(SCALE,{'char'},{'scalartext','nonempty'},...
    'findNearest','SCALE',3);
validateattributes(N,{'numeric'},{'scalar','nonempty','nonnan','finite',...
    'integer','positive','<=',numRowsX},'findNearest','N',4);

switch lower(SCALE)
    case {'linear','l1'}
        distVec = sum(abs(xIn - ones(numRowsX,1)*xD),2);
    case {'least-squares','l2'}
        distVec = sqrt(sum((xIn - ones(numRowsX,1)*xD).^2,2));
    case 'angular'
        distVec = acos(normalizeVector(xIn,2)*normalizeVector(xD,2).');
    otherwise
        error('Invalid specification for SCALE.')
end

[~,indVec] = sort(distVec);
IND = indVec(1:N);
xOut = xIn(IND,:);

end
