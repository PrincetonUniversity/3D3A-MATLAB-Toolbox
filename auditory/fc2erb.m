function erb = fc2erb(fc,n)
%FC2ERB Equivalent rectangular bandwidth (ERB) at a center frequency.
%   ERB = FC2ERB(FC) computes the ERB at center frequency FC, given in Hz.
%
%   ERB = FC2ERB(FC,N) uses the Nth order polynomial approximation given by
%   Moore and Glasberg. Accepts N = 1 or N = 2 only.

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
%   Copyright (c) 2017 Princeton University
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

%   References:
%     [1] Moore and Glasberg (1983) Suggested formulae for calculating
%         auditory-filter bandwidths and excitation patterns.
%     [2] Glasberg and Moore (1990) Derivation of auditory filter shapes
%         from notched-noise data.

if nargin < 2
    n = 1;
end

fc = fc/1000; % convert Hz to kHz

switch n
    case 1
        % The approximation is applicable at moderate sound levels and for
        % values of fc between 0.1 and 10 kHz.
        erb = 24.7*(4.37*fc + 1);
    case 2
        % The approximation is based on the results of a number of
        % published simultaneous masking experiments and is valid from 0.1
        % to 6.5 kHz.
        erb = 6.23*(fc.^2) + 93.39*fc + 28.52;
    otherwise
        error('No polynomial approximation is known for n = %g',n)
end