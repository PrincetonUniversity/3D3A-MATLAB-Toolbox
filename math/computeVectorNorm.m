function AN = computeVectorNorm(A,P,DIM)
%COMPUTEVECTORNORM Compute the norm of a vector.
%   AN = COMPUTEVECTORNORM(A) computes the l2-norm of the vector, A. If A
%   is a matrix, the norm of each column is computed.
%
%   AN = COMPUTEVECTORNORM(A,P) computes the lP-norm. Specify P as inf for
%   the l-infinity norm.
%
%   AN = COMPUTEVECTORNORM(...,DIM) specifies the dimension along which to
%   perform the computation.

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

narginchk(1,3);

if nargin < 3
    DIM = 1;
end

if nargin < 2 || isempty(P)
    P = 2;
end

validateattributes(DIM,{'double'},{'scalar','nonempty','nonnan',...
    'finite','positive','integer','<=',2},'computeVectorNorm','DIM')

if P == inf
    AN = max(abs(A),[],DIM);
else
    AN = sum(abs(A).^P,DIM).^(1/P);
end

end
