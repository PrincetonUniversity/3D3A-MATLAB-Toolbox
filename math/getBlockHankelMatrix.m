function B = getBlockHankelMatrix(A)
%GETBLOCKHANKELMATRIX Generate a block-Hankel matrix.
%   B = GETBLOCKHANKELMATRIX(A) generates a block-Hankel matrix from the 
%   data in A. If the dimensions of A are p x m x N, then the dimensions of
%   B are (p x N) x (m x N).

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

narginchk(1,1);

[p,m,N] = size(A);
B = zeros(p*N,m*N);
BCol1 = zeros(p*N,m);

% Generate entire first column
for ii = 1:N
    p1 = (p*ii)-(p-1);
    BCol1(p1:p1+(p-1),:) = A(:,:,ii);
end

% Populate outputMat using shifted versions of BCol1
for ii = 1:N
    m1 = (m*ii)-(m-1);
    ptd = (ii-1)*p;
    pdt = (p*N)-ptd;
    B(1:pdt,m1:m1+(m-1)) = BCol1((ptd+1):end,:);
end

end
