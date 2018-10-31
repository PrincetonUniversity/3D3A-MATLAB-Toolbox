function P = getPressure(PSI, N, DIM, OPTION)
%GETPRESSURE Inverse discrete Fourier transform for acoustic signals.
%   Computes the pressure due to a potential signal using the convention 
%   shown in GETPOTENTIAL.
%
%   See also GETPOTENTIAL, IFFT.

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
if nargin < 1
    error('No input arguments.');
end

sizeVec = size(PSI);

% No assumption is made about the symmetry of the spectrum by default
if nargin < 4
    OPTION = 'nonsymmetric';
end

% Uses first non-singleton dimension by default
if nargin < 3
    DIM = find(sizeVec ~= 1, 1, 'first');
end

% Uses length along dimension DIM by default
if (nargin < 2) || ((nargin >= 2) && isempty(N))
    N = sizeVec(DIM);
end

P = N*conj(ifft(conj(PSI), N, DIM, OPTION));

end