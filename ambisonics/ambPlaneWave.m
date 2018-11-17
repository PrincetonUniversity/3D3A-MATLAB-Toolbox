function A = ambPlaneWave(L,s0,k,r,t0,ambNorm)
%AMBPLANEWAVE Ambisonics potentials for a plane-wave source.
%   A = AMBPLANEWAVE(L,S0,K,R) computes the ambisonics potentials A, for
%   angular wavenumber K and up to order L, given a plane-wave signal
%   originating from the direction S0 and evaluated at R (both given in
%   Cartesian coordinates).
%
%   A = AMBPLANEWAVE(L,S0,K,R,T0) additionally specifies a the time-of-
%   arrival, T0, at which the plane-wave reaches the origin. By default, no
%   delay is added.
%
%   A = AMBPLANEWAVE(L,S0,K,R,T0,AMBNORM) additionally specifies the
%   ambisonics normalization convention to use. By default, N3D is used.

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

% No delay by default
if nargin < 5 || isempty(t0)
    t0 = 0;
end

% Uses N3D by default
if nargin < 6 || isempty(ambNorm)
    ambNorm = 'N3D';
end

kHat = -normalizeVector(s0); % Propagation direction
psi = exp(1i*k*dot(kHat,r)) .* exp(1i*k*getSoundSpeed()*t0); % Plane-wave potential at r
A = psi * getAmbYMatrix(s0,L,ambNorm).'; % Encoded to ambisonics

end