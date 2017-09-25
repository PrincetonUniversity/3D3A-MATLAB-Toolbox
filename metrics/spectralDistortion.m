function [distSpec,distValF,distValD] = spectralDistortion(hTest,hRef,fS,...
    METHODF,METHODD,SCALE,RANGE)
%SPECTRALDISTORTION Compute absolute distortion spectra and corresponding
%average values for a given input spectrum relative to a given input
%reference spectrum.
%   [distSpec,distValF,distValD] = SPECTRALDISTORTION(hTest,hRef,fS)
%   computes the absolute spectral distortion of hTest relative to hRef and
%   returns the absolute distortion spectrum and distortion values obtained
%   by averaging across rows (frequency average from 20 Hz to 20 kHz 
%   yielding distValF) and also separately across all columns (yielding
%   distValD). Averaging across rows is performed logarithmically by
%   default, whereas across columns, a linear average is calculated.
%   Averaging is performed on the computed distortion spectrum, distSpec,
%   represented in dB. fS is sampling rate in Hz. hTest and hRef may be
%   vectors or 2D matrices of impulse responses (IRs). If matrices, the IRs
%   must be stored as columns.
%
%   [distSpec,distValF,distValD] = SPECTRALDISTORTION(...,METHODF)
%   additionally specifies the method to average distSpec across frequency 
%   to produce distValF. The options are 'log' (logarithmic average, 
%   default), 'lin' (linear average), and 'rms' (root-mean-square).
%
%   [distSpec,distValF,distValD] = SPECTRALDISTORTION(...,METHODD)
%   additionally specifies the method to average distSpec across columns to
%   produce distValD. The options are 'lin' (default), and 'rms'.
%
%   [distSpec,distValF,distValD] = SPECTRALDISTORTION(...,SCALE)
%   additionally specifies the representation of distSpec on which 
%   averaging is performed. The options are 'dB' (magnitude spectrum in dB,
%   default), 'power' (power spectrum), and 'mag' (magnitude spectrum).
%
%   [distSpec,distValF,distValD] = SPECTRALDISTORTION(...,RANGE)
%   additionally specifies the frequency range in Hz over which averaging 
%   should be performed.
%
%   EXAMPLE: Compute spectral distortion of the impulse response h relative
%   to hR, each specified at a sampling rate of fS Hz. Perform rms
%   averaging over frequency ranging from 500 Hz to 15 kHz of the computed 
%   distortion spectrum represented in dB. Also retain the default of
%   averaging linearly across columns.
%   [distSpec,distValF,distValD] = SPECTRALDISTORTION(h,hR,fS,'rms',...
%       'lin','dB',[500,15000]);
%   Note that even though two default values were used, they had to be
%   explicitly specified in order to modify the last value.
%
%   See also DBSPEC, MAGSPEC, COMPUTESPECTRUMAVG, COMPUTEDIRECTIONAVG.

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
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

narginchk(3,7);

if nargin < 7
    RANGE = [20,20000];
end

if nargin < 6
    SCALE = 'dB';
end

if nargin < 5
    METHODD = 'lin';
end

if nargin < 4
    METHODF = 'log';
end

hTest = shiftdim(hTest);
hRef = shiftdim(hRef);
hrirLen = size(hTest,1);
fVec = getFreqVec(fS,hrirLen);

switch lower(SCALE)
    case 'db'
        testDat = dBSpec(hTest);
        refDat = dBSpec(hRef);
        distSpec = abs(refDat-testDat);
        distValF = computeSpectrumAvg(distSpec,fVec,METHODF,RANGE);
        distValD = computeDirectionAvg(distSpec,METHODD);
    case 'power'
        testDat = (magSpec(hTest)).^2;
        refDat = (magSpec(hRef)).^2;
        distSpec = refDat./testDat;
        distValF = computeSpectrumAvg(distSpec,fVec,METHODF,RANGE);
        distValD = computeDirectionAvg(distSpec,METHODD);
    case 'mag'
        testDat = magSpec(hTest);
        refDat = magSpec(hRef);
        distSpec = abs(refDat./testDat);
        distValF = computeSpectrumAvg(distSpec,fVec,METHODF,RANGE);
        distValD = computeDirectionAvg(distSpec,METHODD);
    otherwise
        error('Invalid input for SCALE');
end

end