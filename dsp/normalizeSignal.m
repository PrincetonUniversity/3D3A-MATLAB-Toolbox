function [hNorm,normValdB] = normalizeSignal(h,varargin)
%NORMALIZESIGNAL Normalize a time-domain signal.
%   [HN,NV] = NORMALIZESIGNAL(H) normalizes time-domain signals in H such 
%   that the maximum absolute amplitude of each is 1. If H is a matrix, the
%   signals must be stored as columns and a 'global' normalization is 
%   performed (i.e. the max. absolute amplitude across all signals is set 
%   to 1). The normalized signals in HN are stored as columns. Also
%   returned is the normalization value, NV, in dB.
%
%   ___ = NORMALIZESIGNAL(H,{'f',fVal,'fs',fS}) performs normalization
%   such that the magnitude at frequency fVal (specified in Hz) is 0 dB. 
%   The sampling rate, fS, of the signals in H must be specified in Hz.
%       If fVal is specified as a range of frequencies, fVal = [fL,fU], the 
%       average magnitude in that frequency range is set to 0 dB. 
%       If fVal is set to inf, then the frequency at which the global 
%       maximum magnitude exists is identified and normalization is 
%       performed such that this magnitude is set to 0 dB.
%       If fVal is specified as {inf,[fL,fU]}, then the search for the
%       global maximum magnitude (as described above) is performed within
%       the frequency range [fL,fU].
%
%   ___ = NORMALIZESIGNAL(H,{'f',fVal,'mag',magVal,'fs',fS}) additionally
%   specifies magVal, the magnitude value (in dB) that the normalized 
%   signals must have for any given fVal specification.
%
%   ___ = NORMALIZESIGNAL(H,{...},'local') performs a 'local' normalization
%   (i.e. for a matrix of signals, H, the normalization of each signal in H
%   is performed independently of the other signals in H). If there are no
%   parameters to specify as the second input, use:
%       [HN,NV] = NORMALIZESIGNAL(H,{},'local')

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
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

narginchk(1,3);

% Parse and verify inputs
inputs = parseNormalizeSignalInputs(h,varargin);

% Extract parsed inputs
h = inputs.h;
PARAMS = inputs.PARAMS;
NTYPE = inputs.NTYPE;

flag = 0;
[hLen,numDirs] = size(h);
if isempty(PARAMS)
    refMag = max(abs(h));
    magVal = 0;
else
    indx = find(strcmpi(PARAMS,'f'));
    if isempty(indx)
        error('''f'' must be specified as part of the second input.')
    else
        fVal = PARAMS{indx+1};
    end
    
    indx = find(strcmpi(PARAMS,'mag'));
    if isempty(indx)
        magVal = 0; % Normalize to 0 dB.
    else
        magVal = PARAMS{indx+1};
    end
    
    indx = find(strcmpi(PARAMS,'fs'));
    if isempty(indx)
        error('''fs'' must be specified as part of the second input.')
    else
        fS = PARAMS{indx+1};
    end
    
    % Check specification of optional input values
    validateattributes(fS,{'double'},{'scalar','nonempty','nonnan',...
        'finite','real','positive'},'normalizeSignal','fs');
    validateattributes(magVal,{'double'},{'scalar','nonempty','nonnan',...
        'finite','real'},'normalizeSignal','magVal');
    validateattributes(fVal,{'double','cell'},{'2d','nonempty'},...
        'normalizeSignal','fVal');
    
    fVec = getFreqVec(fS,hLen);
    switch length(fVal)
        case 1
            if fVal == inf
                flag = 1;
                fLIndx = 1;
                fUIndx = ceil((hLen+1)/2);
            elseif (fVal < inf) && (fVal > fS/2)
                error(['fVal must not exceed the Nyquist frequency,',...
                    ' except when specifying fVal as inf.'])
            else
                [~,fLIndx] = min(abs(fVec-fVal));
                fUIndx = fLIndx;
            end
        case 2
            if iscell(fVal)
                flag = 1;
                if fVal{1} == inf
                    [~,fLIndx] = min(abs(fVec-fVal{2}(1)));
                    [~,fUIndx] = min(abs(fVec-fVal{2}(2)));
                else
                    error('If fVal is a cell array, fVal{1} must be inf.')
                end
            else
                if fVal(1) > fVal(2)
                    error('If fVal = [fL,fU], fL must be <= fU.')
                else
                    [~,fLIndx] = min(abs(fVec-fVal(1)));
                    [~,fUIndx] = min(abs(fVec-fVal(2)));
                end
            end
        otherwise
            error('fVal can only either be a scalar or 2-element array.')
    end
    
    hMagdB = getMagSpecdB(h);
    hMagdB_crop = hMagdB(fLIndx:fUIndx,:);
    if flag == 1
        [~,fIndx] = max(hMagdB_crop);
        refMag = zeros(1,numDirs);
        for ii = 1:numDirs
            refMag(ii) = db2mag(mean(hMagdB_crop(fIndx(ii),:)));
        end
    else
        refMag = db2mag(mean(hMagdB_crop,1));
    end
end

switch lower(NTYPE)
    case 'global'
        normVal = db2mag(magVal)/max(refMag);
    case 'local'
        normVal = repmat(db2mag(magVal)./refMag,hLen,1);
    otherwise
        error('Invalid specification for third input.')
end

hNorm = h.*normVal;
normValdB = mag2db(normVal);

end

function inputs = parseNormalizeSignalInputs(h,opts)
%PARSENORMALIZESIGNALINPUTS Parse and verify inputs to the normalizeSignal
%function.

p = inputParser;

% Required inputs
addRequired(p,'h',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite'},'normalizeSignal','h',1));
h = shiftdim(h); % If h is a vector, force to a column.

% Optional inputs
addOptional(p,'PARAMS',{},@(x)validateattributes(x,{'cell'},...
    {'vector'},'normalizeSignal','parameters for second input',2));
addOptional(p,'NTYPE','global',@(x)validateattributes(x,{'char'},...
    {'scalartext'},'normalizeSignal','parameter for third input',3));

p.CaseSensitive = false;
p.FunctionName = 'normalizeSignal';

parse(p,h,opts{:});

inputs = p.Results;

end
