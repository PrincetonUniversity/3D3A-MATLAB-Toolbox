function PSI = getPotential(P, N, DIM)
%GETPOTENTIAL Discrete Fourier transform for acoustic signals.
%   PSI = GETPOTENTIAL(P) computes the potential due to a pressure signal
%   using the following convention:
% 
%   For a length N pressure signal P, the potential is a length N vector
%   PSI, with elements
%                       N
%       PSI(k) = (1/N) sum  P (n)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
%                      n=1
%   The inverse transform (computed by getPressure) is given by
%                       N
%        P (n) =       sum PSI(k)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N.
%                      k=1
%
%   PSI = GETPOTENTIAL(P,N) computes the N-point FFT.
%
%   PSI = GETPOTENTIAL(P,N,DIM) computes the FFT along dimension DIM.
%
%   PSI = GETPOTENTIAL(P,[],DIM) computes the FFT of length SIZE(P,DIM).
%
%   See also GETPRESSURE, FFT.

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

% Needs at least 1 input argument
narginchk(1,3);

sizeVec = size(P);

% Uses first non-singleton dimension by default
if nargin < 3
    DIM = find(sizeVec ~= 1, 1, 'first');
end

% Uses length along dimension DIM by default
if (nargin < 2) || isempty(N)
    N = sizeVec(DIM);
end

PSI = (1/N)*conj(fft(conj(P), N, DIM));

end