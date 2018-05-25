function hInv = inverseFilter(h,reg,varargin)
%INVERSEFILTER Compute an inverse filter with optional regularization.
%   hInv = INVERSEFILTER(h) returns the inverse of h without any
%   regularization. h may be a vector or 2D matrix with the data to be
%   inverted stored as columns. hInv has the same dimensions as h. This
%   command may also be specified as: hInv = INVERSEFILTER(h,'none');
%
%   hInv = INVERSEFILTER(h,'piecewise') returns the inverse of h with
%   frequency-dependent regularization using a piecewise regularization 
%   profile. Options for specifying a custom regularization profile are 
%   listed below. If a custom regularization profile is not provided, no
%   regularization is performed.
%
%   hInv = INVERSEFILTER(h,'piecewise','range',range) specifies the regularization
%   frequency range and regularization value for each "piece" of the
%   piecewise regularization function as a 3x1 row vector [w1,w2,eps]. 
%   Multiple "pieces" are specified as separate rows. eps is the 
%   regularization parameter in the range [w1,w2]. w1 and w2 are normalized
%   frequencies with 0 corresponding to DC and 1 the Nyquist frequency. The
%   default is [0,1,0] (i.e. no regularization).
%       Example: range = [0,0.1,0.001;0.2,0.8,0;0.9,1,0.001] is a "shelf" 
%       profile with both low and high shelves. Transitions between 0.1 and
%       0.2, and 0.8 and 0.9 are linear.
%
%   hInv = INVERSEFILTER(h,'gardner','sampleRate',fS) returns the inverse 
%   of h with frequency-dependent regularization described by Gardner
%   and Martin [1]. Sampling rate, fS, must also be specified, in Hz, as 
%   a value-pair. By default, regularization is performed from 1 octave to
%   7 octaves above the smallest, non-zero frequency. The dynamic range of
%   the inverse filter is also limited to 24 dB, by default.
%
%   hInv = INVERSEFILTER(h,'gardner','sampleRate',fS,Name1,Value1,...) 
%   specifies optional comma-separated pairs of Name,Value arguments, 
%   where Name is the argument name and Value is the corresponding value. 
%   Name must appear inside single quotes (' '). You can specify two name 
%   and value pair arguments in any order as Name1,Value1,Name2,Value2.
%   Valid Name,Value arguments are as follows:
%
%   'avgRange'      1x2 row vector specifying frequency range over which 
%                   direct inversion is performed. Both frequencies must be
%                   specified in Hz. Default: [1 octave above lowest freq.,
%                   7 octaves above lowest freq.].
%
%   'dynRange'      Inverse filter dynamic range in dB (default: 24).
%                   Minimum value is 1 dB. Maximum value is 150 dB.
%
%   EXAMPLE: Design an inverse filter from h using frequency-dependent
%   regularization described by Gardner and Martin [1] with default
%   settings except for inverse filter dynamic range, which is now 20 dB.
%       hInv = INVERSEFILTER(h,'gardner','sampleRate',fS,'dynRange',20);

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

%   References:
%       [1] Gardner and Martin (1994) - HRTF Measurements of a KEMAR Dummy 
%       Head Microphone.pdf

narginchk(1,8);

[hLen,numCh] = size(h); % assumes input signals stored as columns
if isvector(h)
    h = shiftdim(h); % forces h to a column if it's a row vector.
elseif ~ismatrix(h)
    error('Signal to be inverted must be a vector or matrix.');
end

if nargin < 2
    reg = 'none'; % no regularization
end

if strcmpi(reg,'piecewise')
    ind = find(strcmpi(varargin,'range'),1);
    if isempty(ind)
        range = [0,1,0];
    else
        range = varargin{ind+1};
        if ~ismatrix(range)
            error('range must be a 3xn matrix with n >= 1.')
        end
    end
elseif strcmpi(reg,'gardner')
    params = {'sampleRate','avgRange','dynRange'};
    indVec = cell(length(params),1);
    for ii = 1:length(params)
        indVec{ii} = find(strcmpi(varargin,params{ii}),1);
    end
    if isempty(indVec{1})
        error('sampleRate must be specified for reg -> gardner')
    else
        sampleRate = varargin{indVec{1}+1};
    end
    if isempty(indVec{2})
        fVec = getFreqVec(sampleRate,hLen);
        lowFreq = 2*fVec(2); % 1 octave above lowest non-zero freq.
        highFreq = (2^7)*fVec(2); % 7 octaves above lowest non-zero freq.
        highFreq = clip(highFreq,[lowFreq,20000]);
        avgRange = [lowFreq,highFreq];
    else
        avgRange = varargin{indVec{2}+1};
        if length(avgRange) ~= 2
            error('avgRange must be a row vector with 2 entries.')
        end
    end
    if isempty(indVec{3})
        dynRange = 24;
    else
        dynRange = varargin{indVec{3}+1};
    end
end

switch lower(reg)
    case {'piecewise'}
        H = fft(h,[],1);
        err = zeros(ceil((hLen+1)/2),1);
        lfIndVec = round(range(:,1)*(ceil((hLen+1)/2)-1))+1;
        hfIndVec = round(range(:,2)*(ceil((hLen+1)/2)-1))+1;
        numPieces = size(range,1);
        for ii = 1:numPieces-1
            err(lfIndVec(ii):hfIndVec(ii),1) = range(ii,3);
            if lfIndVec(ii+1) > hfIndVec(ii)
                m = (range(ii+1,3)-range(ii,3))/...
                    (lfIndVec(ii+1)-hfIndVec(ii));
                err(hfIndVec(ii)+1:lfIndVec(ii+1)-1,1) = range(ii,3)+...
                    m*((hfIndVec(ii)+1:lfIndVec(ii+1)-1).'-hfIndVec(ii));
            end
        end
        err(lfIndVec(numPieces):hfIndVec(numPieces),1) = range(numPieces,3);
        fullErr = repmat([err;flipud(err(2:(end-1),1))],1,numCh);
        hInv = ifft(conj(H)./((conj(H).*H)+fullErr),[],1,'symmetric');
    case {'gardner'}
        nyqFreq = sampleRate/2;
        [~,nyqInd] = min(abs(fVec-nyqFreq));
        minLowFreq = fVec(1);
        maxHighFreq = fVec(nyqInd);
        avgRange = clip(avgRange,[minLowFreq,maxHighFreq]);
        dynRange = clip(dynRange,[1,150]);
        H = fft(h,[],1);
        HdBMag = mag2db(abs(H));
        avgdBMag = repmat(logmean(HdBMag,fVec,avgRange),hLen,1);
        dynHalfRange = dynRange/2;
        Hcomp = HdBMag - avgdBMag;
        Hcomp(Hcomp > dynHalfRange) = dynHalfRange;
        Hcomp(Hcomp < -dynHalfRange) = -dynHalfRange;
        Hcomp = Hcomp + avgdBMag;
        HInv = db2mag(-Hcomp).*exp(-1i*angle(H));
        hInv = ifft(HInv,[],1,'symmetric');
    otherwise
        hInv = ifft(1./fft(h),[],1,'symmetric');
end

hInv = reshape(hInv,[hLen,numCh]);

end

