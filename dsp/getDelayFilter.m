function h = getDelayFilter(FFTLen, Fs, t, ADJUST)
%GETDELAYFILTER Sinc-filter delay impulse response.
%   H = GETDELAYFILTER(N,FS,T) returns the N-point impulse response H given
%   at sampling rate FS for a time-delay of T seconds.
%
%   H = GETDELAYFILTER(N,FS,T,ADJUST) optionally adjusts the time-delay to
%   to an integer number of samples if ADJUST evaluates to true.

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

if nargin < 4 || isempty(ADJUST)
    ADJUST = false;
end

if ADJUST
    d = mod(round(t*Fs),FFTLen);
    h = zeros(FFTLen,1);
    h(d+1) = 1;
else
    specLen = 1 + FFTLen/2;
    H = exp(-1i*2.0*pi*t*Fs/FFTLen*(0:specLen).');
    h = ifft(H,FFTLen,1,'symmetric');
end

end