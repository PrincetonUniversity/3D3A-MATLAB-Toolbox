function p = a2p(a,Fs,r,ambNorm)
%A2P Convert ambisonics signals to acoustic pressure.
%   P = A2P(A,FS,R) computes the acoustic pressure P at a position R (given
%   in Cartesian coordinates) due to a spherical Fourier-Bessel series
%   expansion with time-domain signals A at sampling rate FS. If R is empty
%   or omitted, the origin [0 0 0] is used.
%   
%   A should be a K-by-(L+1)^2 matrix where K is the signal length and L is
%   the maximum (truncation) order of the expansion. In this case, P will
%   be a K-by-1 column vector.
%
%   P = A2P(A,FS,R,AMBNORM) additionally specifies the ambisonics
%   normalization convention used. By default, N3D is assumed.
%
%   See also A2PSI, GETPOTENTIAL, GETPRESSURE.

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

narginchk(2,4);

% Renders at the origin by default
if nargin < 3 || isempty(r)
    r = [0 0 0];
end

% Uses N3D normalization by default
if nargin < 4 || isempty(ambNorm)
    ambNorm = 'N3D';
end

FFTLen = size(a,1);
A = getPotential(a,FFTLen,1);
k = f2k(getFreqVec(Fs,FFTLen));
PSI = A2psi(A,k,r,ambNorm);
p = getPressure(PSI,FFTLen,1,'symmetric');

end