function [ILD, fc, avgILD] = estimateILD(hL,hR,Fs,B,FRANGE)
%ESTIMATEILD Interaural level difference spectra.
%   ILD = ESTIMATEILD(HL,HR) estimates the interaural level difference
%   (ILD) spectra for the binaural (or head-related) impulse responses HL
%   and HR. HL and HR must be the same size, and ILD will be a matrix also 
%   of the same size.
%
%   [ILD,FC] = ESTIMATEILD(HL,HR,FS) returns the frequency vector FC.
%   corresponding to the rows of ILD, when the sampling rate FS is
%   specified. If FS is omitted or empty, FC will be a normalized frequency
%   with 1.0 corresponding to the Nyquist frequency and 2.0 corresponding
%   to the sampling rate.
%
%   [ILD,FC,AVG] = ESTIMATEILD(HL,HR,FS) additionally returns AVG, the
%   logmean average of ILD over all non-zero frequencies.
%
%   [ILD,FC,AVG] = ESTIMATEILD(HL,HR,FS,B,[FMIN,FMAX]) if either the
%   bandwidth B (in octaves) or the frequency range [FMIN,FMAX] (in Hz) is
%   specified, then ILD will be returned as band-averaged spectra, with FC
%   corresponding to the center frequencies of each band and AVG computed
%   by taking a linear mean over all bands.
%
%   Alternatively, B may be specified as 'erb' to use ERB-spaced auditory
%   bands (implemented as a gammatone filterbank).
%
%   See also COMPUTEBANDAVG, GETGAMMATONEFILTERS.

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

narginchk(2,5);

FFTLen = size(hL,1);

bandAvgFlag = false;
if nargin < 4 || isempty(B)
    B = 'erb';
else
    bandAvgFlag = true;
end

if nargin < 5 || isempty(FRANGE)
    FRANGE = [20 20000];
else
    bandAvgFlag = true;
end

if nargin >= 3 && ~isempty(Fs)
    f = getFreqVec(Fs,FFTLen);
else
    f = getFreqVec(2,FFTLen); % normalized frequency: Fs = 2, Nyquist = 1
    if bandAvgFlag
        warning('Cannot compute band-averages without a sampling rate.');
        bandAvgFlag = false;
    end
end

rawILD = mag2db(abs(fft(hL,FFTLen,1)./fft(hR,FFTLen,1)));
if bandAvgFlag
    [ILD, fc] = computeBandAvg(rawILD,f,B,FRANGE,Fs);
    avgILD = mean(ILD,1);
else
    ILD = rawILD;
    fc = f;
    avgILD = logmean(rawILD,f);
end

end
