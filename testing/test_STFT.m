% Testing forward and inverse STFT

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

xLen = 1000;
winLen = 256;
noverlap = 128;
nfft = 512;
padflag = true;

windowNames = {'rectwin','hann','hamming','blackman'};
windows = cell(length(windowNames),1);
z = cell(length(windowNames),1);

figure()
hold all
for ww = 1:length(windowNames)
    % Generate input signal
    x = [0; randn(xLen-1,1)];
    % Note: for window functions that start with 0, the first sample cannot
    % be reconstructed.
    
    % Compute window function
    eval(['windows{ww} = ' windowNames{ww} '(winLen+1);']);
    windows{ww} = windows{ww}(1:winLen);
    
    % Compute STFT
    Y = getForwardSTFT(x, windows{ww}, noverlap, nfft, padflag);
    
    % Compute inverse STFT
    z{ww} = getInverseSTFT(Y, windows{ww}, noverlap, nfft, padflag);
    
    % Plot discrepancy
    plot(x-z{ww}(1:xLen))
end
legend(windowNames)