function Y = getForwardSTFT(x, window, noverlap, nfft)
%GETFORWARDSTFT Spectrogram using short-time Fourier transform (STFT).
%   Y = GETFORWARDSTFT(X,WINDOW,NOVERLAP) returns Y, the STFT of a signal
%   X, using the specified WINDOW vector and overlapping NOVERLAP samples.
%
%   Y = GETFORWARDSTFT(X,WINDOW,NOVERLAP,NFFT) computes NFFT-length FFTs at
%   each time frame. If unspecified, NFFT = LENGTH(WINDOW).
%
%   See also SPECTROGRAM, GETINVERSESTFT.

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

winLen = length(window);
if nargin < 4
    nfft = winLen;
end

x = shiftdim(x);
xLen = length(x);
hopLen = winLen - noverlap;
numPartitions = len2part(xLen + winLen, noverlap, hopLen);
xPadLen = part2len(numPartitions, winLen, hopLen);

Y = spectrogram([zeros(hopLen,1); x; zeros(xPadLen - (xLen + hopLen),1)], window, noverlap, nfft, 'twosided'); % NFFT x numTimeFrames

end

function nparts = len2part(len, novlp, hop)
nparts = ceil((len - novlp) / hop);
end

function len = part2len(nparts, wlen, hop)
len = wlen + hop * (nparts - 1);
end