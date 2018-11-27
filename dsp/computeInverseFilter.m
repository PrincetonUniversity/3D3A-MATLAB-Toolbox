function hInv = computeInverseFilter(h,varargin)
%COMPUTEINVERSEFILTER Compute an inverse filter.
%   hInv = COMPUTEINVERSEFILTER(h) returns the inverse of h without any
%   regularization. h may be a vector, or matrix with the data to be
%   inverted stored as columns. If h is a vector, hInv is a column vector 
%   with the same length as h. If h is a matrix, hInv has the same 
%   dimensions as h. This command may also be specified as: 
%       hInv = COMPUTEINVERSEFILTER(h,'direct');
%
%   ___ = COMPUTEINVERSEFILTER(h,'piecewise') returns the inverse of h with
%   frequency-dependent regularization using a piecewise regularization 
%   profile, with regularization performed outside the range [w1,w2], where 
%   w1 = 0.1 and w2 = 0.9, are normalized frequencies (normalized with 
%   respect to the Nyquist frequency). The amount of regularization, eps, 
%   applied is 0.001.
%
%   ___ = COMPUTEINVERSEFILTER(h,'piecewise',PARAMS) specifies parameters 
%   for regularization. PARAMS specifies the regularization frequency range 
%   and regularization value for each "piece" of the piecewise 
%   regularization function as a 3x1 row vector [w1,w2,eps]. Multiple 
%   "pieces" are specified as separate rows with eps being the amount of 
%   regularization applied in the range [w1,w2]. w1 and w2 are normalized 
%   frequencies with 0 corresponding to DC and 1 the Nyquist frequency.
%   EXAMPLES:
%       1. [0,1,0] - no regularization.
%       2. [0,0.1,0.001; 0.2,0.8,0; 0.9,1,0.001] - "shelf" profile with 
%       both low and high shelves. Here, regularization is performed below 
%       a normalized frequency of 0.1 and above a normalized frequency of 
%       0.8. Transitions between 0.1 and 0.2, and 0.8 and 0.9 are linear. 
%       No regularization is performed between 0.2 and 0.8 since eps is set 
%       to 0 in this range.
%
%   ___ = COMPUTEINVERSEFILTER(h,'gardner1994') returns the inverse of h 
%   with frequency-dependent regularization described by Gardner and Martin 
%   [1]. Regularization is performed from 1 octave to 7 octaves above the
%   smallest, non-zero frequency, subject to a maximum of the Nyquist
%   frequency. The dynamic range of the inverse filter is also limited to
%   24 dB.
%
%   ___ = COMPUTEINVERSEFILTER(h,'gardner1994',PARAMS) specifies parameters 
%   for regularization. PARAMS must be a cell array in which the following 
%   may be specified as Name-Value pairs.
%       1. 'avgRange',[w1,w2] - vector specifying frequency range over
%           which direct inversion is performed. As before, w1 and w2 are 
%           normalized frequencies.
%       2. 'dynRange',DR - DR specifies the dynamic range, in dB, of the  
%           inverse filter subject to: 1 <= DR <= 150.
%
%   See also GETPIECEWISEPROFILE.

%   =======================================================================
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
%       [1] Gardner and Martin (1994) - HRTF Measurements of a KEMAR Dummy 
%       Head Microphone.

narginchk(1,3);

% Parse and verify inputs
[inputs,extra] = parseComputeInverseFilterInputs(h,varargin);

% Extract parsed inputs
h = inputs.h;
TYPE = inputs.TYPE;
PARAMS = inputs.PARAMS;

h = shiftdim(h); % If h is a vector, force it to be a column vector.
[hLen,numChs] = size(h);
switch lower(TYPE)
    case 'direct'
        hInv = ifft(1./fft(h),'symmetric');
    case 'piecewise'
        H = fft(h);
        nyqIndx = ceil((hLen+1)/2);
        pProfileHalf = getPiecewiseProfile(PARAMS,nyqIndx);
        pProfile = repmat([pProfileHalf;...
            flipud(pProfileHalf(2:(end-1),1))],1,numChs);
        hInv = ifft(conj(H)./((conj(H).*H)+pProfile),'symmetric');
    case 'gardner1994'
        % Extract avgRange from PARAMS
        indx = find(strcmpi(PARAMS,'avgRange'),1);
        if isempty(indx)
            avgRange = extra{1};
        else
            avgRange = PARAMS{indx+1};
        end
        % Verify avgRange
        validateattributes(avgRange,{'double'},{'vector','nonempty',...
            'nonnan','finite','real','nonnegative','numel',2,'<=',1},...
            'computeInverseFilter','avgRange')
        
        % Extract dynRange from PARAMS
        indx = find(strcmpi(PARAMS,'dynRange'),1);
        if isempty(indx)
            DR = extra{2};
        else
            DR = PARAMS{indx+1};
        end
        % Verify DR
        validateattributes(DR,{'double'},{'scalar','nonempty','nonnan',...
            'finite','real','positive','<=',150},'computeInverseFilter',...
            'DR specification for dynRange')
        
        wVec = getFreqVec(2,hLen);
        w1 = min(avgRange);
        w2 = max(avgRange);
        H = getMagSpecdB(h);
        Hp = getPhaseSpec(h);
        avgMagdB = repmat(logmean(H,wVec,[w1,w2]),hLen,1);
        halfDR = DR/2;
        H_compressed = H - avgMagdB;
        H_compressed(H_compressed > halfDR) = halfDR;
        H_compressed(H_compressed < -halfDR) = -halfDR;
        H_compressed = H_compressed + avgMagdB;
        HInv = db2mag(-H_compressed).*exp(-1i*Hp);
        hInv = ifft(HInv,'symmetric');
    otherwise
        error('Unrecognized second input.')
end

end

function [inputs,extra] = parseComputeInverseFilterInputs(h,opts)
%PARSECOMPUTEINVERSEFILTERINPUTS Parse and verify inputs for the 
%computeInverseFilter function.

p = inputParser;

% Required inputs
addRequired(p,'h',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite','real'},'computeInverseFilter','h',1));
h = shiftdim(h); % If h is a vector, force it to be a column vector.

% Optional inputs
addOptional(p,'TYPE','direct',@(x)validateattributes(x,{'char'},...
    {'scalartext','nonempty'},'computeInverseFilter','type of inverse'));
extra = {}; % No additional parameters need to be returned.
if ~isempty(opts)
    switch lower(opts{1})
        case 'piecewise'
            defaultVal = [0,0.1,0.001; 0.2,0.8,0; 0.9,1,0.001];
            addOptional(p,'PARAMS',defaultVal,@(x)validateattributes(x,...
                {'double'},{'2d','nonempty','nonnan','finite','real',...
                'ncols',3},'computeInverseFilter',['PARAMS for option:',...
                ' ''piecewise''']));
        case 'gardner1994'
            hLen = size(h,1);
            w1 = 2*(2/hLen); % 1 octave above smallest, non-zero freq.
            w2 = (2^7)*(2/hLen); % 7 octaves above smallest, non-zero freq.
            if w1 >= 1
                w1 = 0;
                warning(['Length of h is very short. Setting w1 = 0 ',...
                    'for avgRange.'])
            end
            if w2 >= 1
                w2 = 1;
                warning(['Length of h is very short. Setting w2 = 1 ',...
                    'for avgRange.'])
            end
            DR = 24; % Default dynamic range in dB.
            defaultVal = {'avgRange',[w1,w2],'dynRange',DR};
            addOptional(p,'PARAMS',defaultVal,@(x)validateattributes(x,...
                {'cell'},{'nonempty','nrows',1},'computeInverseFilter',...
                'PARAMS for option: ''gardner1994'''));
            extra = {[w1,w2],DR}; % Default values returned for performing  
                                  % additional checks.
        otherwise
            addOptional(p,'PARAMS',[]); % Unused, so no validation needed.
    end
else
    addOptional(p,'PARAMS',[]); % Unused, so no validation needed.
end

p.CaseSensitive = false;
p.FunctionName = 'computeInverseFilter';

parse(p,h,opts{:});

inputs = p.Results;

end
