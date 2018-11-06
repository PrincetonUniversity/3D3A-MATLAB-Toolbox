function a = mu2a(mu,r,L,w,ambNorm)
%MU2A Convert plane-wave signals to ambisonics signals.
%   A = MU2A(MU,R,L) computes ambisonics signals A, up to order L, given
%   plane-wave signals MU or source directions R. MU should be a matrix of
%   size K-by-Q, where Q is the number of source directions; A will be a
%   matrix of size K-by-(L+1)^2. R should be a matrix of size Q-by-3.
%
%   A = MU2A(MU,R,L,W) additionally specifies the quadrature weights W. By
%   default, equal weights are used, with W = 4*PI/Q. Note that SUM(W)
%   should equal approximately 4*PI.
%
%   A = MU2A(MU,R,L,W,AMBNORM) additionally specifies the ambisonics
%   normalization convention to use. By default, N3D is used.
%
%   See also A2MU, GETAMBYMATRIX.

%   ==============================================================================
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

narginchk(3,5);

nCols = size(r,1);

if nargin < 4 || isempty(w)
    w = ones(nCols,1)*(4*pi/nCols);
end

if nargin < 5 || isempty(ambNorm)
    ambNorm = 'N3D';
end

Y = getAmbYMatrix(r, L, ambNorm);
a = mu*diag(w)*(Y.');

end