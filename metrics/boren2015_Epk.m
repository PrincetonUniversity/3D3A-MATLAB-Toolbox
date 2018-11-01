function [Epk, En] = boren2015_Epk(h, h0, Fs, FRANGE, SCALE, BW, THRESH)
%BOREN2015_EPK Boren's peak and notch errors.
%   [EPK,EN] = BOREN2015_EPK(H,H0,FS) computes peak and notch errors, EPK
%   and EN, respectively, given an input signal H and reference signal H0,
%   given at sampling rate FS.
%
%   [EPK,EN] = BOREN2015_EPK(H,[],FS) uses a flat reference spectrum.
%
%   [EPK,EN] = BOREN2015_EPK(H,H0,FS,FRANGE) computes errors over the
%   specified frequency range FRANGE (by default, [50 21000] is used).
%
%   [EPK,EN] = BOREN2015_EPK(H,H0,FS,FRANGE,SCALE) uses the specified
%   frequency scale, which can be either 'lin' (or 'hz') for linear, 'oct'
%   for octaves (default), and 'erb' for equivalent rectangular bandwidth.
%
%   [EPK,EN] = BOREN2015_EPK(H,H0,FS,FRANGE,SCALE,BW) specifies that peaks
%   must be at least BW apart on the specified SCALE.
%       For 'lin' (or 'hz'), BW = 100 Hz by default.
%       For 'oct', BW = 1/3 octave by default.
%       For 'erb', BW = 1 ERB by default.
%
%   [EPK,EN] = BOREN2015_EPK(H,H0,FS,FRANGE,SCALE,BW,THRESH) uses the
%   specified threshold (in dB) for peak detection. By default, THRESH = 1
%   dB.
%
%   See also FINDPEAKS, FRACTIONALOCTAVESMOOTH.

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
%     [1] Boren et al. (2015) Coloration metrics for headphone equalization.

FFTLen = length(h);
freqVec = getFreqVec(Fs, FFTLen);

if nargin < 4
    FRANGE = [50 21000]; % 50 Hz to 21 kHz by default
end

if nargin < 5 || isempty(SCALE)
    SCALE = 'oct'; % octave scale by default
end

if nargin < 6
    switch lower(SCALE)
        case {'hz','lin'}
            BW = 100; % 100 Hz bandwidths by default
        case 'oct'
            BW = 1/3; % 1/3-octave bandwidths by default
        case 'erb'
            BW = 1; % 1 ERB bandwidths by default
    end
end

if nargin < 7
    THRESH = 1; % 1 dB threshold by default
end

if isempty(h0) % reference = flat
    H = fft(h, FFTLen);
else
    H = fft(h, FFTLen)./fft(h0, FFTLen);
end
Hf = fractionalOctaveSmooth(H, 48);
Hc = fractionalOctaveSmooth(H,  1);
Hd = mag2db(abs(Hf./Hc));

[~, ind1] = findNearest(freqVec, FRANGE(1));
[~, ind2] = findNearest(freqVec, FRANGE(2));
indx = ind1:ind2;

switch lower(SCALE)
    case {'hz','lin'}
        scaleFreq = freqVec(indx);
    case 'oct'
        scaleFreq = log2(freqVec(indx)/FRANGE(1));
    case 'erb'
        if exist('freqtoaud','file') == 2
            scaleFreq = freqtoaud(freqVec(indx),'erb'); % Needs LTFAT
        else
            error('freqtoaud from LTFAT not found.');
        end
end

[E1, ~] = findpeaks( Hd(indx), scaleFreq, 'MinPeakHeight', THRESH, 'MinPeakDistance', BW/2);
[E2, ~] = findpeaks(-Hd(indx), scaleFreq, 'MinPeakHeight', THRESH, 'MinPeakDistance', BW/2);
% [E1, ~] = findDistinctSpectralPeaks( Hd(indx), freqVec(indx), dbThresh, bwOctave);
% [E2, ~] = findDistinctSpectralPeaks(-Hd(indx), freqVec(indx), dbThresh, bwOctave);

nbands = range(scaleFreq)/BW;
Epk = sum(E1)/nbands;
En = sum(E2)/nbands;

end

function [Yp, Xp] = findDistinctSpectralPeaks(Y, X, Ymin, bwOctave)

Ys = 0;
Xs = 0;
if any(Y >= Ymin)
    [Y1, X1] = findpeaks(Y,X,'MinPeakHeight',Ymin);
    [Ys, Is] = sort(Y1,'descend');
    Xs = X1(Is);
    
    bw = 2^(bwOctave/2);
    ii = 1;
    while ii <= length(Ys)
        temp = (Ys < Ys(ii)) & (Xs < Xs(ii)*bw) & (Xs > Xs(ii)/bw);
        Ys(temp) = [];
        Xs(temp) = [];
        ii = ii + 1;
    end
end
Yp = Ys;
Xp = Xs;

end