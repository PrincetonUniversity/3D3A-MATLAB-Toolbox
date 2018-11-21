function Y = getAmbYMatrix(r, L, ambNorm)
%GETAMBYMATRIX Matrix of real-valued spherical harmonics.
%   Y = GETAMBYMATRIX(R,L,AMBNORM) returns the (L+1)^2-by-SIZE(R,1) matrix
%   of real-valued spherical harmonics, evaluated in the directions R up to
%   order L and with normalization convention AMBNORM. By default, N3D
%   normalization is used.
%
%   See also AMBSPHERICALHARMONICY.

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

narginchk(2,3);

% Uses N3D normalization by default
if nargin < 3 || isempty(ambNorm)
    ambNorm = 'N3D';
end

nRows = size(r,1);
N = (L + 1)^2;

Y = zeros(N,nRows);
for n = 0:(N-1)
    [l, m] = getAmbOrder(n);
    Y(n+1,:) = ambSphericalHarmonicY(l, m, r, ambNorm);
end

end