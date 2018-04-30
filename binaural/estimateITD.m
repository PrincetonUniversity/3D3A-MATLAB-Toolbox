function ITD = estimateITD(hL,hR,Fs,method,varargin)
%ESTIMATEITD Estimate the interaural time difference (ITD), in seconds.
%   ITD = ESTIMATEITD(hL,hR,Fs) estimates ITD for input HRIRs hL and hR at 
%   sampling rate Fs specified in Hz using thresholding with a threshold of
%   0.2 (20% of the absolute maximum of each IR). Negative ITDs correspond
%   to sound sources on the left. hL and hR may be vectors or matrices but
%   must have the same dimensions. If matrices, the IRs must be stored in
%   columns.
%
%   ITD = ESTIMATEITD(...,method) additionally specifies the method of 
%   estimating ITD. The options are:
%   1. 'group delay' or 'group' estimates ITD accurate to 1 sample at 
%   the specified sampling rate by computing the difference in the group
%   delays of hL and hR, and averaging the result over nominal frequencies 
%   of 0 to 1500 Hz. To modify the range over which averaging is performed,
%   see the additional options below.
%   2. 'cross-correlation' or 'xcorr' estimates ITD accurate to 1 sample
%   at the specified sampling rate by determining the first sample position 
%   at which the cross-correlation between hL and hR reaches an absolute
%   maximum.
%   3. 'phase delay' or 'phase' estimates ITD with fractional-sample 
%   accuracy by determing the difference between the unwrapped phase at the
%   right ear and that at the left ear, each averaged over 0 to 1500 Hz
%   nominal frequencies. To modify the range over which averaging is 
%   performed, see the additional options below.
%   4. 'thresholding' or 'thresh' estimates ITD accurate to 1 sample at the
%   specified sampling rate by thresholding with a threshold of 0.2 (20% of
%   the absolute maximum of each IR). This is the default.
%
%   ITD = ESTIMATEITD(...,Name1,Value1,...) specifies optional 
%   comma-separated pairs of Name,Value arguments, where Name is the 
%   argument name and Value is the corresponding value. Name must appear 
%   inside single quotes (' '). You can specify the following:
%   
%   'filter'        Specifications for a Butterworth filter to filter hL
%                   and hR prior to ITD estimation. The specifications must
%                   be a 1x3 cell array specifying {filter order, cut-off
%                   frequency in Hz, filter type} for a Butterworth
%                   filter. For order and type specifications, see the
%                   documentation for 'butter'. The three parameters must
%                   be specified in the order shown above. 'zerophase' is 
%                   an optional fourth parameter that may be specified to
%                   design a zero-phase filter.
%
%   'resample'      Resample hL and hR prior to estimating ITD. This must
%                   be a scalar > 0 such that a value in the range (0,1)
%                   corresponds to downsampling, and a value > 1
%                   corresponds to upsampling. This corresponds to
%                   specifying parameter 'P' in the 'resample' function
%                   assuming 'Q' is set to 1.
%
%   'range'         Frequency range, in Hz, over which averaging is
%                   performed for the 'group delay' and 'phase delay' 
%                   estimation methods. This must be specified as a 1x2 or
%                   a 2x1 vector with the first value specifying the lower
%                   bound.
%
%   'threshold'     Threshold fraction used for the thresholding approach
%                   to estimate ITD. This must be a scalar between 0 and 1,
%                   with 1 corresponding to the maximum absolute value of a
%                   given IR, and 0 corresponding to the first non-zero
%                   sample.
%
%   'linearfit'     [No Value] Specifying this parameter first generates a
%                   linear fit to the unwrapped phase response over the 
%                   frequency range specified in 'range' (or the default of
%                   0 to 1500 Hz) prior to estimating ITD using the 'phase 
%                   delay' approach.

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
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

if nargin < 3
    error('Not enough input arguments.')
end
if ~isequal(size(hL),size(hR))
    error('Input HRIRs must be the same size.')
end

% Perform optional filtering first
indx = findInCell(varargin,'filter');
if indx
    specCell = varargin{indx+1};
    n = specCell{1,1}; % filter order; see 'butter' help
    cutoff = specCell{1,2}; % cutoff frequency in Hz
    ftype = specCell{1,3}; % filter type; see 'butter' help
    Wn = 2*pi*cutoff/Fs;
    [b,a] = butter(n,Wn,ftype);
    if length(specCell) == 4 && strcmpi(specCell{1,4},'zerophase')
        hL = filtfilt(b,a,hL);
        hR = filtfilt(b,a,hR);
    else
        hL = filter(b,a,hL);
        hR = filter(b,a,hR);
    end
end

% Perform resampling next ('upsample' included for backwards-compatibility)
indx = findInCell(varargin,'resample')+findInCell(varargin,'upsample');
if indx
    Fs = varargin{indx+1}*Fs;
    hL = resample(hL,varargin{indx+1},1);
    hR = resample(hR,varargin{indx+1},1);
end

if isvector(hL)
    hL = shiftdim(hL);
    hR = shiftdim(hR);
    FFTLen = length(hL);
    numIRs = 1;
else
    [FFTLen, numIRs] = size(hL);
end
halfLen = ceil((FFTLen+1)/2);

freqVec = getFreqVec(Fs,FFTLen);
indx = findInCell(varargin,'range');
if indx
    if numel(varargin{indx+1})==1
        fU = varargin{indx+1};
        fL = 0;
    else
        fL = varargin{indx+1}(1);
        fU = varargin{indx+1}(2);
    end
else
    fL = 0;
    fU = 1500;
end
[~,nL] = findNearest(freqVec,fL);
[~,nU] = findNearest(freqVec,fU);
[nU,~] = min([nU;halfLen]); % Prevents choosing value above Nyquist freq.
avgRange = nL:nU;

indx = findInCell(varargin,'threshold');
if indx
    thp = varargin{indx+1};
else
    thp = 0.2;
end

switch lower(method)
    case {'group delay','group'}
        dL = real(fft(diag(0:FFTLen-1)*hL,FFTLen,1)./fft(hL,FFTLen,1));
        dR = real(fft(diag(0:FFTLen-1)*hR,FFTLen,1)./fft(hR,FFTLen,1));
        d = mean(dL(avgRange,:)-dR(avgRange,:),1,'omitnan');
    case {'cross-correlation','xcorr'}
        d = zeros(1,numIRs);
        for ii = 1:numIRs
            [xc,lagVec] = xcorr(hL(:,ii),hR(:,ii));
            [~,lagIndex] = max(abs(xc));
            d(ii) = lagVec(lagIndex);
        end
    case {'phase delay','phase'}
        pL = unwrap(angle(fft(hL,FFTLen,1)),[],1);
        pR = unwrap(angle(fft(hR,FFTLen,1)),[],1);
        indx = findInCell(varargin,'linearfit');
        if indx
            for ii = 1:numIRs
                fitL = polyfit(freqVec(avgRange),pL(avgRange,ii),1);
                fitR = polyfit(freqVec(avgRange),pR(avgRange,ii),1);
                pL(:,ii) = polyval(fitL,freqVec);
                pR(:,ii) = polyval(fitR,freqVec);
            end
        end
        d = Fs*mean(diag(1./(2*pi*freqVec(avgRange)))*(pR(avgRange,:)-...
            pL(avgRange,:)),1,'omitnan');
    otherwise % Thresholding
        dL = thresholdIRs(hL,thp);
        dR = thresholdIRs(hR,thp);
        d = dL-dR;
end

ITD = d/Fs;

end
