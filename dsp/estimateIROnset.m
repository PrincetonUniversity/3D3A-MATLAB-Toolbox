function onsetMat = estimateIROnset(inputIR,varargin)
%ESTIMATEIRONSET Estimate onset of an impulse response.
%   onsetMat = ESTIMATEIRONSET(inputIR) estimates the onset of an impulse
%   response in inputIR using thresholding with a 20% threshold. inputIR
%   can be a vector or matrix. If it's a matrix, the IRs should be stored
%   as columns. The returned onset is specified in samples.
% 
%   ___ = ESTIMATEIRONSET(...,METHOD) optionally specifies the method to
%   use to estimate onsets. METHOD must be specified as a cell array with 
%   the format {methodName,methodParameters}. The following may be
%   specified for METHOD:
%       1. {'threshold',thp} where thp is the thresholding percentage to
%       use. thp can take values in the range [0,100].
%       2. {'grpdelay',fS} where fS is the sampling rate in Hz. This
%       estimates onset as the average group delay between 0 Hz and the
%       Nyquist frequency for the specified sampling rate.
%           2.i. {'grpdelay',fS,[fL,fU]} optionally averages the group
%           delay in the range [fL,fU], where fL and fU are in Hz.
%       3. {'phase',fS} where fS is the sampling rate in Hz. This
%       estimates onset as the average slope of the unwrapped phase 
%       response between 0 Hz and the Nyquist frequency for the specified
%       sampling rate.
%           3.i. {'phase',fS,[fL,fU]}, optionally averages the slope in the
%           range [fL,fU], where fL and fU are in Hz.
%           3.ii. {'phase',fS,[fL,fU],'linearfit'} optionally applies a
%           linear fit to the unwrapped phase response before computing the
%           average slope. The fit is also applied over [fL,fU].
%       4. {'xcorr'} estimates onset as the sample value corresponding to
%       the max. absolute value of the cross-correlation spectrum of the
%       inputIR and its minimum-phase version.
%
%   See also IRGRPDELAY.

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
inputs = parseESTIMATEIRONSETInputs(inputIR,varargin);

% Extract parsed inputs
inputIR = inputs.inputIR;
METHOD = inputs.METHOD;

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
            'thp for METHOD{1} = ''threshold''')
        if thp <= 1
            warning(['Entered value of thp is <= 1. Note that thp is ',...
                'specified as a percentage ranging from 0 to 100.'])
        end
        onsetMat = thresholdIRs(inputIR,thp/100);
    case 'grpdelay'
        if length(METHOD) > 1
            fS = METHOD{2};
            validateattributes(fS,{'double'},{'scalar','real','finite',...
            'nonnan','positive'},'estimateIROnset',['fS for ',...
            'METHOD{1} = ''grpdelay'''])
        else
            error('Sampling rate in Hz must be specified.');
        end
        if length(METHOD) < 3
            avgRange = [0,fS/2];
        else
            avgRange = METHOD{3};
        end
        validateattributes(avgRange,{'double'},{'vector','real',...
            'finite','nonnan','nonnegative','numel',2,'<=',fS/2},...
            'estimateIROnset','avgRange for METHOD{1} = ''grpdelay''')
        [~,onsetMat] = irGrpDelay(inputIR,fS,avgRange);
    case 'phase'
        if length(METHOD) > 1
            fS = METHOD{2};
            validateattributes(fS,{'double'},{'scalar','real','finite',...
            'nonnan','positive'},'estimateIROnset',['fS for METHOD{1}',...
            ' = ''phase'''])
        else
            error('Sampling rate in Hz must be specified.');
        end
        fVec = getFreqVec(fS,irLen);
        nyqFreq = fS/2;
        nyqIndx = find(fVec == nyqFreq,1);
        if length(METHOD) < 3
            avgRange = [0,nyqFreq];
        else
            avgRange = METHOD{3};
        end
        validateattributes(avgRange,{'double'},{'vector','real',...
            'finite','nonnan','nonnegative','numel',2,'<=',nyqFreq},...
            'estimateIROnset','avgRange for METHOD{1} = ''phase''')
        [~,fLIndx] = min(abs(fVec-avgRange(1)));
        [~,fUIndx] = min(abs(fVec-avgRange(2)));
        fUIndx = min([fUIndx,nyqIndx]);
        if length(METHOD) > 3
            validateattributes(METHOD{4},{'char'},{'scalartext'},...
            'estimateIROnset','METHOD{4} for METHOD{1} = ''phase''')
            if strcmpi(METHOD{4},'linearfit')
                fitFlag = 1;
            else
                error(['Unrecognized optional input for METHOD{1}',...
                    ' = ''phase''.']);
            end
        else
            fitFlag = 0;
        end
        phaseSpec = unwrap(angle(fft(inputIR)));
        if fitFlag
            for ii = 1:numIRs
                fittingLine = polyfit(fVec(fLIndx:fUIndx),...
                    phaseSpec(fLIndx:fUIndx,ii),1);
                phaseSpec(:,ii) = polyval(fittingLine,fVec);
            end
        end
        onsetMat = -fS*mean(diag(1./(2*pi*fVec(fLIndx:fUIndx)))*...
            phaseSpec(fLIndx:fUIndx,:),'omitnan');
    case 'xcorr'
        onsetMat = zeros(1,numIRs);
        minPhaseIR = makeMinPhaseIR(inputIR,'hilb');
        for ii = 1:numIRs
            [xc,lagVec] = xcorr(inputIR(:,ii),minPhaseIR(:,ii));
            [~,lagIndex] = max(abs(xc));
            onsetMat(ii) = lagVec(lagIndex);
        end
    otherwise
        error('Invalid method specification');
end

end

function inputs = parseESTIMATEIRONSETInputs(inputIR,opts)
%PARSEESTIMATEIRONSETINPUTS Parse and verify inputs for the estimateIROnset 
%function.

p = inputParser;

% If inputIR is a vector, force it to be a column vector.
if isvector(inputIR)
    inputIR = shiftdim(inputIR);
end

% Required inputs
addRequired(p,'inputIR',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite','real'},'estimateIROnset','inputIR',1));

% Optional inputs
addOptional(p,'METHOD',{'threshold',20},@(x)validateattributes(x,...
    {'cell'},{'nonempty','size',[1,NaN]},'estimateIROnset','METHOD'));

p.CaseSensitive = false;
p.FunctionName = 'estimateIROnset';

parse(p,inputIR,opts{:});

inputs = p.Results;

end
