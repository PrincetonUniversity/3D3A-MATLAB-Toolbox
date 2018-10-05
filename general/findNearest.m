function [x1,IND] = findNearest(xVec,x0,SCALE,N)
%FINDNEAREST Find the row of a matrix nearest to a given row vector.
%   X1 = FINDNEAREST(X,X0) returns X1, the row of the matrix X which is
%   closest in value to the row vector X0.
%
%   [X1,I1] = FINDNEAREST(X,X0) returns the position I1 of X1 in X.
%
%   [X1,I1] = FINDNEAREST(X,X0,SCALE) measures distance between elements
%   using the specified SCALE (see options below).
%
%   [X1,I1] = FINDNEAREST(X,X0,SCALE,N) returns N nearest elements.
%   
%   OPTIONAL INPUTS:
%   1. SCALE (default = 'linear') - specifies the distance metric to use. 
%   The options are 'linear' or 'l1', 'least-squares' or 'l2', and
%   'angular'.
%   2. N (default = 1) - specifies the number of elements to return.

%   ==============================================================================
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
%   Copyright (c) 2018 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   ==============================================================================

narginchk(2,4);

% Re-orient arrays if needed
if isrow(xVec) && ~isequal(size(xVec), size(x0))
    xVec = shiftdim(xVec);
end

% Use linear distance metric by default
if nargin < 3 || isempty(SCALE)
    SCALE = 'linear';
end

% Find 1 nearest point by default
if nargin < 4 || isempty(N)
	N = 1;
end

[xLen, ~] = size(xVec);
switch lower(SCALE)
    case {'linear', 'l1'}
        distVec = sum(abs(xVec - ones(xLen,1)*x0),2);
    case {'least-squares', 'l2'}
        distVec = sqrt(sum((xVec - ones(xLen,1)*x0).^2,2));
    case 'angular'
        distVec = acos(normalizeVector(xVec,2)*normalizeVector(x0,2).');
end
[~, indVec] = sort(distVec);
IND = indVec(1:N);
x1 = xVec(IND,:);

end