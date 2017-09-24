function avgVal = computeSpectrumAvg(inputSpectra,fVec,METHOD,RANGE)
%COMPUTESPECTRUMAVG Average an input spectrum over frequency.
%   avgVal = COMPUTESPECTRUMAVG(inputSpectra,fVec) computes the logarithmic 
%   average, over frequency, of the spectra in inputSpectra. fVec is a 
%   vector of frequencies specified in Hz. The frequencies in fVec 
%   correspond to the frequencies at which spectral data in inputSpectra 
%   are specified. inputSpectra may be a vector or 2D matrix. If it's a 
%   matrix, individual spectra must be stored as columns. avgVal is either 
%   a scalar or vector depending on whether inputSpectra is a vector or 2D 
%   matrix, respectively.
%
%   avgVal = COMPUTESPECTRUMAVG(...,METHOD) additionally specifies the 
%   method used to average the data in inputSpectra across frequency to 
%   produce avgVal. The options are 'log' (default), 'lin' (linear average)
%   , and 'rms' (root-mean-square).
%
%   avgVal = COMPUTESPECTRUMAVG(...,METHOD,RANGE) additionally specifies 
%   the frequency range, in Hz, over which averaging should be performed.
%
%   See also LOGMEAN, MEAN, RMS.

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

narginchk(2,4);

if nargin < 4
    RANGE = [20,20000];
end

if nargin < 3
    METHOD = 'log';
end

fVec = shiftdim(fVec);

[~,lInd] = findNearest(fVec,RANGE(1));
[~,hInd] = findNearest(fVec,RANGE(2));

switch lower(METHOD)
    case 'log'
        avgVal = logmean(inputSpectra,fVec,RANGE);
    case 'lin'
        avgVal = mean(inputSpectra(lInd:hInd,:));
    case 'rms'
        avgVal = rms(inputSpectra(lInd:hInd,:));
    otherwise
        error('Invalid input for METHOD.');
end

end

