function mu = a2mu(a,r,ambNorm,pinvFlag,w)
%A2MU Convert ambisonics signals to plane-wave signals.
%   MU = A2MU(A,R) computes the plane-wave signals MU (also called the
%   signature function) for source directions R given ambisonics signals A.
%   A should be a matrix of size K-by-N, where N is the number of
%   ambisonics channels; MU will be a matrix of size K-by-Q, where Q is the
%   number of source directions. R should be a matrix of size Q-by-3.
%
%   MU = A2MU(A,R,AMBNORM) additionally specifies the ambisonics
%   normalization convention to use. By default, N3D is used.
%
%   MU = A2MU(A,R,AMBNORM,PINVFLAG) additionally specifies the conversion
%   method to use, where if PINVFLAG evaluates to true, a least-squares
%   inversion is taken to compute MU. By default, PINVFLAG is false.
%
%   MU = A2MU(A,R,AMBNORM,PINVFLAG,W) additionally specifies quadrature
%   weights W for the grid of directions R. By default equal weights are
%   used, with W = 4*PI/Q. Note that SUM(W) should equal approximately 4*PI.
%
%   See also MU2A, GETAMBYMATRIX.

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

narginchk(2,5);

N = size(a,2);
Q = size(r,1);

if nargin < 3 || isempty(ambNorm)
    ambNorm = 'N3D';
end

if nargin < 4 || isempty(pinvFlag)
    pinvFlag = false;
end

if nargin < 5 || isempty(w)
    w = ones(Q,1)*(4*pi/Q);
    if pinvFlag
        warning(['Using pinv method but no quadrature weights were '...
            'specified; using equal weights.']);
    end
end

Y = getAmbYMatrix(r, sqrt(N) - 1, ambNorm);
if pinvFlag
    mu = (a/(Y.'))*diag(1./w);
else % beamforming by default
    F = zeros(N,1);
    for n = 0:(N-1)
        [l, ~] = getAmbOrder(n);
        F(n+1) = ambNormSquared(l, ambNorm);
    end
    mu = a*diag(1./F)*Y;
end

end