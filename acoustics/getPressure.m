function P = getPressure(varargin)
%GETPRESSURE Inverse discrete Fourier transform for acoustic signals.
%   P = GETPRESSURE(PSI) computes the pressure P due to a potential signal 
%   PSI.
%
%   P = GETPRESSURE(...,'symmetric') returns a real-valued signal P by
%   forcing conjugate-symmetry of PSI.
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
if ischar(varargin{end})
    maxNargin = 4;
    OPTION = varargin{end};
else
    maxNargin = 3;
    OPTION = 'nonsymmetric'; % No assumption is made about the symmetry of the spectrum by default
end
narginchk(1,maxNargin);

PSI = varargin{1};
sizeVec = size(PSI);
switch maxNargin - nargin
    case 0
        DIM = varargin{3};
        N = varargin{2};
    case 1
        DIM = find(sizeVec ~= 1, 1, 'first'); % Uses first non-singleton dimension by default
        N = varargin{2};
    case 2
        DIM = find(sizeVec ~= 1, 1, 'first');
        N = sizeVec(DIM); % Uses length along dimension DIM by default
end

P = N*conj(ifft(conj(PSI), N, DIM, OPTION));

end