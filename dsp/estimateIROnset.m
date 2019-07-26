function onsetMat = estimateIROnset(inputIR,varargin)
%ESTIMATEIRONSET Estimate onset of an impulse response (IR).
%   O = ESTIMATEIRONSET(X) estimates the onset of impulse response(s), X,
%   using thresholding with a 20% relative threshold.
%       If X is a vector, O will be a scalar.
%       If X is a matrix containing N IRs, the IRs should be stored as
%       columns. O will then be a row vector of length N.
%   The returned onset, O, is specified in samples and is accurate to 1 
%   sample.
%
%   O = ESTIMATEIRONSET(X,METHOD) optionally specifies the method to use
%   when estimating onsets. METHOD must be specified as a cell array with
%   the format {Name,Parameters}. The following may be specified:
%       1. {'threshold',THP} where THP is the thresholding percentage to
%       use. THP can take values in the range [0,100]. The returned sample 
%       onset, O, is accurate to 1 sample.
%       2. {'grpdelay',Fs} where Fs is the sampling rate in Hz. This
%       estimates onset as the average group delay computed between 0 and
%       Fs/2 Hz.
%           2.1. {'grpdelay',Fs,[FL,FU]} optionally averages the group
%           delay in the range [FL,FU], where FL and FU are in Hz.
%       The returned sample onset, O, is accurate to fractions of 1 sample.
%       3. {'phase',Fs} where Fs is the sampling rate in Hz. This estimates
%       onset as the average slope of the unwrapped phase response between
%       0 and Fs/2 Hz.
%           3.1. {'phase',Fs,[FL,FU]}, optionally averages the slope in the
%           range [FL,FU], where FL and FU are in Hz.
%           3.2. {'phase',Fs,[FL,FU],'linearfit'} optionally applies a
%           linear fit to the unwrapped phase response before computing the
%           slope. The fit is also applied over the range [FL,FU].
%       The returned sample onset, O, is accurate to fractions of 1 sample.
%       4. {'mpxc'} estimates onset as the sample value corresponding to
%       the max. absolute value of the cross-correlation spectrum of X and
%       its minimum-phase version computed using the MAKEMINPHASEIR
%       function. The returned sample onset, O, is accurate to 1 sample.
%
%   Needs: Signal Processing Toolbox.
%
%   See also THRESHOLDIRS, GETGRPDELAY.

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

narginchk(1,2);

% Parse and verify inputs
inputs = parseEstimateIROnsetInputs(inputIR,varargin);

% Extract parsed inputs
inputIR = inputs.inputIR;
METHOD = inputs.METHOD;

% Perform additional input verification
validateattributes(METHOD{1},{'char'},{'scalartext','nonempty'},...
    'estimateIROnset','the first element of METHOD',2)

[irLen,numIRs] = size(inputIR);
switch lower(METHOD{1})
    case 'threshold'
        if length(METHOD) < 2
            thp = 20;
        else
            thp = METHOD{2};
        end
        validateattributes(thp,{'double'},{'scalar','real','finite',...
            'nonnan','nonnegative','<=',100},'estimateIROnset',...
            'the THP value when METHOD Name is ''threshold''',2)
        onsetMat = thresholdIRs(inputIR,thp/100);
    case 'grpdelay'
        if length(METHOD) > 1
            Fs = METHOD{2};
            validateattributes(Fs,{'double'},{'scalar','real','finite',...
                'nonnan','positive'},'estimateIROnset',['the value ',...
                'for Fs ','when METHOD Name is ''grpdelay'''],2)
        else
            error(['Sampling rate in Hz must be specified as an input',...
                ' when METHOD Name is ''grpdelay''.']);
        end
        if length(METHOD) < 3
            avgRange = [0,Fs/2];
        else
            avgRange = METHOD{3};
        end
        validateattributes(avgRange,{'double'},{'vector','real',...
            'finite','nonnan','nonnegative','numel',2,'<=',Fs/2},...
            'estimateIROnset',['the values for [FL,FU] when METHOD ',...
            'Name is ''grpdelay'''],2)
        [~,onsetMat] = getGrpDelay(inputIR,Fs,avgRange);
    case 'phase'
        if length(METHOD) > 1
            Fs = METHOD{2};
            validateattributes(Fs,{'double'},{'scalar','real','finite',...
                'nonnan','positive'},'estimateIROnset',['the value',...
                ' for Fs when METHOD Name is ''phase'''],2)
        else
            error(['Sampling rate in Hz must be specified as an input',...
                ' when METHOD Name is ''phase''.']);
        end
        fVec = getFreqVec(Fs,irLen);
        nyqFreq = Fs/2;
        if length(METHOD) < 3
            avgRange = [0,nyqFreq];
        else
            avgRange = METHOD{3};
        end
        validateattributes(avgRange,{'double'},{'vector','real',...
            'finite','nonnan','nonnegative','numel',2,'<=',nyqFreq},...
            'estimateIROnset',['the values for [FL,FU] when METHOD ',...
            'Name is ''phase'''],2)
        [~,fLIndx] = min(abs(fVec-avgRange(1)));
        [~,fUIndx] = min(abs(fVec-avgRange(2)));
        if length(METHOD) > 3
            validateattributes(METHOD{4},{'char'},{'scalartext'},...
                'estimateIROnset',['the ''linearfit'' option when',...
                ' METHOD Name is ''phase'''],2)
            if strcmpi(METHOD{4},'linearfit')
                fitFlag = 1;
            else
                error(['Unrecognized parameter when METHOD Name is ',...
                    'specified as ''phase''.']);
            end
        else
            fitFlag = 0;
        end
        phaseSpec = getPhaseSpec(inputIR,'unwrap');
        if fitFlag
            for ii = 1:numIRs
                fittingLine = polyfit(fVec(fLIndx:fUIndx),...
                    phaseSpec(fLIndx:fUIndx,ii),1);
                phaseSpec(:,ii) = polyval(fittingLine,fVec);
            end
        end
        onsetMat = -Fs*mean(diag(1./(2*pi*fVec(fLIndx:fUIndx)))*...
            phaseSpec(fLIndx:fUIndx,:),'omitnan');
    case 'mpxc'
        onsetMat = zeros(1,numIRs);
        minPhaseIR = makeMinPhaseIR(inputIR,'hilb');
        for ii = 1:numIRs
            [xc,lagVec] = xcorr(inputIR(:,ii),minPhaseIR(:,ii));
            [~,lagIndex] = max(abs(xc));
            onsetMat(ii) = lagVec(lagIndex);
        end
    otherwise
        error('Invalid METHOD Name specification.');
end

end

function inputs = parseEstimateIROnsetInputs(inputIR,opts)
%PARSEESTIMATEIRONSETINPUTS Parse and verify inputs for the estimateIROnset
%function.

p = inputParser;

% Required inputs
addRequired(p,'inputIR',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite','real'},'estimateIROnset','X',1));

% If inputIR is a vector, force it to be a column vector.
inputIR = shiftdim(inputIR);

% Optional inputs
addOptional(p,'METHOD',{'threshold',20},@(x)validateattributes(x,...
    {'cell'},{'nonempty','size',[1,NaN]},'estimateIROnset','METHOD',2));

p.CaseSensitive = false;
p.FunctionName = 'estimateIROnset';

parse(p,inputIR,opts{:});

inputs = p.Results;

end
