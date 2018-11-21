function [CSdB, fc, FT, CB, A] = kates1984(x, Fs, fc)
%KATES1984 Kates' central spectrum coloration model.
%   CS = KATES1984(X,FS) computes Kates' central spectrum CS (in dB) for a
%   mono signal X given at sampling rate FS.
%
%   [CS,FC,FT,CB,A] = KATES1984(X,FS) additionally returns center
%   frequencies FC, Fourier-transform coefficients FT, critical-band
%   energies CB, and the reference level A.
%
%   [...] = KATES1984(X,FS,FC) uses specified center frequencies FC.

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

%   References:
%     [1] Kates (1984) A Perceptual Criterion for Loudspeaker Evaluation.
%     [2] Kates (1985) A central spectrum model for the perception of 
%         coloration in filtered Gaussian noise.

if nargin < 3
    fc = shiftdim(erbspacebw(20,20000));
else
    fc = shiftdim(fc);
end
numfc = length(fc);

R = xcorr(x); % Autocorrelation function

tau = 70e-3; % 70 ms
nu = Fs*tau;
xLen = length(x);
n = shiftdim((-xLen+1):(xLen-1));
W = R.*exp(-abs(n)/nu); % Weighted autocorrelation function

N = length(R);
t = n/Fs;
FT = zeros(numfc,1);
for kk = 1:numfc
    FT(kk) = sqrt(dot(W,cos(2*pi*fc(kk)*t))); % Fourier-transform coefficient
end
% Confirmed -- matches Kates (1985) Fig. 3

PS = shiftdim(abs(fft(circshift(W,xLen)))); % Power spectrum
f = shiftdim(getFreqVec(Fs,N)).';
B = fc2erb(fc,2); % ERB from Moore and Glasberg (1983)
b = B/(2*((2^(1/3) - 1)^(1/2)));
CB = zeros(numfc,1);
for kk = 1:numfc
    HPS = abs(1 + ((f-fc(kk))/b(kk)).^2).^(-3); % Critical-band filter bank power spectra
%     CB(kk) = sqrt(dot(HPS,PS)/(N*B(kk)/Fs)); % Critical-band energy,
%     "normalized by the filter bandwidths" (Kates, 1984)
    CB(kk) = sqrt(dot(HPS,PS)/sum(HPS)); % Critical-band energy
end
% Confirmed -- matches Kates (1985) Fig. 4

w = 2*pi*fc;
fh = 70; ph = -2*pi*fh; HPF1 = 1i*w./(1i*w - ph); % High-pass filter at 70 Hz
fl = 3000; pl = -2*pi*fl; LPF1 = pl./(1i*w - pl); % Low-pass filter at 3 kHz
ASC = abs(HPF1.*LPF1); % Auditory sensitivity curve
X_LE = sum((CB.^2).*ASC); % Loudness estimate
A = sqrt(X_LE/sum(ASC)); % Reference level

fl = 500; p = -2*pi*fl;
LP = p./(1i*w - p); % Low-pass filter at 500 Hz

a = 0.5;
CS2 = a*(abs(CB).^2 - A^2) + (1-a)*(abs(LP).^2).*(abs(FT).^2 - A^2) + A^2;
CSdB = mag2db(sqrt(CS2)); % Central spectrum
% Confirmed -- matches Kates (1985) Figs. 5 and 6 (when plots are normalized)

% Fs = 100000;
% fc = shiftdim(2*logspace(1,4,256));
% x = zeros(2048,1);
% x(1) = 1;
% x1 = x; x1(101) = 0.5; x1 = x1/sum(abs(x1));
% [C1, ~] = KATES1984(x1,Fs,fc);
% x2 = x; x2(101) = -0.5; x2 = x2/sum(abs(x2));
% [C2, ~] = KATES1984(x2,Fs,fc);
% plot(fc,C1); hold all; plot(fc,C2);
% set(gca,'XScale','log')

end