function [Sc_dB, A0, SD] = wittek2007(b, b0, Fs, fc)
%WITTEK2007 Wittek's spectral alterations coloration model.
%   SC = WITTEK2007(B,B0,FS) computes the spectral alterations (in dB)
%   between binaural signals B and reference binaural signals B0, specified
%   at sampling rate FS.
%
%   SC = WITTEK2007(B,[],FS) uses a flat reference spectrum by default.
%
%   [SC,A0,SD] = WITTEK2007(B,B0,FS) additionally returns Wittek's A0
%   measure and spectral deviation SD.
%
%   [...] = WITTEK2007(B,B0,FS,FC) uses specified center frequencies FC.
%
%   See also GETPATTERSONFILTERS.

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

%   References:
%     [1] Salomons (1995) Coloration and Binaural Decoloration of Sound due
%         to Reflections.
%     [2] Wittek et al. (2007) On the sound colour properties of wavefield
%         synthesis and stereo.

if nargin < 4 || isempty(fc)
    fc = erbspacebw(20,20000);
end

if isempty(b0)
    refFlag = false;
    IRLen = size(b,1);
    Sc0_sq = 1;
else
    refFlag = true;
    IRLen = max(size(b,1),size(b0,1));
end

FFTLen = 2^nextpow2(IRLen);
specLen = 1 + FFTLen/2;

% Compute binaural power spectra
B_sq = abs(fft(b,FFTLen,1)).^2; % FFTLen-by-2
if refFlag % Compute reference binaural power spectra
    B0_sq = abs(fft(b0,FFTLen,1)).^2; % FFTLen-by-2
end

% Compute critical band filter bank
h = getPattersonFilters(fc, Fs, FFTLen);
H = fft(h,FFTLen,1);
Cfc = abs(H(1:specLen,:)); % specLen-by-numfc

% Compute internal spectrum
Sc_sq = diag(1./sum(Cfc,1))*(Cfc.'*B_sq(1:specLen,:)); % Salomons, Eq. 5.12
Sc_sq = mean(Sc_sq,2); % numfc-by-1
if refFlag % Compute reference internal spectrum
    Sc0_sq = diag(1./sum(Cfc,1))*(Cfc.'*B0_sq(1:specLen,:)); % Salomons, Eq. 5.12
    Sc0_sq = mean(Sc0_sq,2); % numfc-by-1
end

% Compute spectral alterations compared to reference
Sc_dB = mag2db(sqrt(Sc_sq./Sc0_sq));
A0 = wittek2007_A0(Sc_dB);
SD = wittek2007_SD(Sc_dB, fc);

end

function A0 = wittek2007_A0(Sc_dB)
%wittek2007_A0 Wittek's A0 measure.

A0 = max(Sc_dB) - min(Sc_dB);

end

function SD = wittek2007_SD(Sc_dB, fc)
%wittek2007_SD Wittek's spectral deviation.

f = logspace(log10(min(fc)), log10(max(fc)));
Sc_dB = interp1(fc,Sc_dB,f);

SD = std(Sc_dB);

end