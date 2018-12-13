function [Y, fc] = computeBandAvg(Q,F,B,FRANGE,Fs)
%COMPUTEBANDAVG Averages of a function in specified frequency bands.
%   Y = COMPUTEBANDAVG(Q,F,B,[FMIN,FMAX]) computes the log-weighted
%   averages of Q, whose values are specified as a function of F, for bands
%   of width B (in octaves) within FMIN and FMAX. F should have uniformly
%   (linearly) spaced values. Q must have the same number of rows as F.
%
%   Alternatively, B may be specified as 'erb' to use ERB-spaced auditory
%   bands (implemented as a gammatone filterbank), in which case the
%   sampling rate Fs must also be specified:
%   Y = COMPUTEBANDAVG(Q,F,'erb',[FMIN,FMAX],Fs)
%   and F must be given as subset of a standard frequency vector.
%
%   See also LOGMEAN, GETGAMMATONEFILTERS, GETFREQVEC.

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

narginchk(4,5);

if isnumeric(B) % Fractional-octave bands
    numfc = 1 + round(log2(FRANGE(2)/FRANGE(1))/B);
    fc = zeros(numfc,1);
    
    Y = zeros(numfc,size(Q,2));
    for ii = 1:numfc
        fc(ii) = FRANGE(1)*(2^((ii-1)*B));
        fl = fc(ii)/(2^(B/2));
        fu = fc(ii)*(2^(B/2));
        Y(ii,:) = logmean(Q, F, [fl, fu]);
    end
elseif strcmpi(B,'erb') % ERB-spaced auditory bands
    if exist('erbspacebw','file') == 2
        fc = erbspacebw(FRANGE(1),FRANGE(2));
        FFTLen = round(Fs/mean(diff(F)));
        
        f = getFreqVec(Fs,FFTLen);
        h = getGammatoneFilters(fc, Fs, FFTLen);
        Hmag = getMagSpec(h,1); % FFTLen-by-numfc
        
        % Resample to frequencies given in F
        indx = ismember(f,F);
        HmagF = Hmag(indx,:); % length(F)-by-numfc
        
        Y = (diag(1./sum(HmagF,1))*(HmagF.'))*Q;
    else
        warning('erbspacebw from LTFAT not found, using 1/3-octave bands instead.');
        [Y, fc] = computeBandAverage(Q,F,1/3,FRANGE,Fs);
    end
end

end