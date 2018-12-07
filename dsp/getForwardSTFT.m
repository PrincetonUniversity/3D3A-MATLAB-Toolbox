function [Y,T] = getForwardSTFT(x, window, noverlap, nfft, padFlag)
%GETFORWARDSTFT Spectrogram using short-time Fourier transform (STFT).
%   Y = GETFORWARDSTFT(X,WINDOW,NOVERLAP) returns Y, the STFT of a signal
%   X, using the specified WINDOW vector and overlapping NOVERLAP samples.
%
%   Y = GETFORWARDSTFT(X,WINDOW,NOVERLAP,NFFT) computes NFFT-length FFTs at
%   each time frame. If unspecified, NFFT = LENGTH(WINDOW).
%
%   [Y,T] = GETFORWARDSTFT(...) additionally returns the time position (in
%   samples) of the center of each partition.
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
if nargin < 4 || isempty(nfft)
    nfft = winLen;
end
if nargin < 5 || isempty(padFlag)
    padFlag = false;
end

x = shiftdim(x);
xLen = length(x);
hopLen = winLen - noverlap;
if padFlag
    numPartitions = STFT_len2part(xLen, winLen, noverlap, padFlag);
    xPadLen = STFT_part2len(numPartitions, winLen, noverlap, false); % return padded length
    x = [zeros(hopLen,1); x; zeros(xPadLen - (xLen + hopLen),1)];
end

Y = spectrogram(x, window, noverlap, nfft, 'twosided'); % NFFT x numPartitions
T = ((winLen/2):hopLen:(length(x)-winLen/2)).' - double(padFlag)*hopLen;

end