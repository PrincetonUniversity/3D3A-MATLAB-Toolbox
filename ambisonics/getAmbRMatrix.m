function V = getAmbRMatrix(k, L, r, ambNorm)
%GETAMBRMATRIX Matrix of regular Helmholtz solutions.
%   V = GETAMBRMATRIX(K,L,R,AMBNORM) returns the LENGTH(K)-by-(L+1)^2
%   matrix V of regular solutions to the 3D Helmholtz equation, evaluated
%   at a single position R (given in Cartesian coordinates) for angular
%   wavenumbers K and up to order L, with normalization convention AMBNORM.
%   By default, N3D normalization is used.
%
%   See also AMBHELMHOLTZSOLUTIONR.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Joseph G. Tylka <josephgt@princeton.edu>
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

narginchk(3,4);

% Uses N3D normalization by default
if nargin < 4 || isempty(ambNorm)
    ambNorm = 'N3D';
end

nRows = length(k);
N = (L + 1)^2;

V = zeros(nRows,N);
for n = 0:(N-1)
    [l, m] = getAmbOrder(n);
    V(:,n+1) = ambHelmholtzSolutionR(l, m, shiftdim(k), r, ambNorm);
end

end