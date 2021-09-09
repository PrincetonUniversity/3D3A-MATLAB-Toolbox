function varargout = estimateIROnset(inputIR,varargin)
%ESTIMATEIRONSET Onset of an impulse response (IR).
%   O = ESTIMATEIRONSET(X) estimates the onset of impulse response(s), X,
%   using thresholding with a 20% relative threshold.
%       If X is a vector, O will be a scalar.
%       If X is a matrix containing N IRs, the IRs should be stored as
%       columns. O will then be a row vector of length N.
%   The returned onset, O, is specified in samples and is accurate to 1 
%   sample unless the IRs are resampled (see option below).
%
%   O = ESTIMATEIRONSET(X,METHOD) optionally specifies the method to use
%   when estimating onsets. METHOD must be specified as a cell array with
%   the format {Name,Parameters}. The following may be specified:
%       1. {'threshold',THP} where THP is the thresholding percentage to
%       use. THP can take values in the range [0,100]. The returned sample 
%       onset, O, is accurate to 1 sample unless the IRs are resampled (see
%       option below).
%
%       2. {'grpdelay',Fs} where Fs is the sampling rate in Hz. This
%       estimates onset as the average group delay computed between 0 and
%       Fs/2 Hz.
%           2.1. {'grpdelay',Fs,[FL,FU]} optionally averages the group
%           delay in the range [FL,FU], where FL and FU are in Hz with
%           0 <= FL <= FU <= Fs/2.
%       The returned sample onset, O, is accurate to fractions of 1 sample.
%
%       3. {'phase',Fs} where Fs is the sampling rate in Hz. This estimates
%       onset as the average slope of the unwrapped phase response between
%       0 and Fs/2 Hz.
%           3.1. {'phase',Fs,[FL,FU]}, optionally averages the slope in the
%           range [FL,FU], where FL and FU are in Hz with 
%           0 <= FL <= FU <= Fs/2.
%           3.2. {'phase',Fs,[FL,FU],'linearfit'} optionally applies a
%           linear fit to the unwrapped phase response before computing the
%           slope. The fit is applied over the range [FL,FU].
%       The returned sample onset, O, is accurate to fractions of 1 sample.
%
%       4. {'mpxc'} estimates onset as the sample value corresponding to
%       the max. absolute value of the cross-correlation spectrum of X and
%       its minimum-phase version computed using the 'makeMinPhaseIR'
%       function. The returned sample onset, O, is accurate to fractions of 
%       1 sample.
%
%       5. {'mpxc2',Fs} estimates onset as the sample value corresponding   
%       to the max. absolute value of the cross-correlation spectrum of a 
%       windowed version of X and its minimum-phase version computed using 
%       the 'makeMinPhaseIR' function. The returned sample onset, O, is 
%       accurate to fractions of 1 sample. Sampling rate, Fs, must be 
%       specified in Hz.
%
%       6. {'multithresh',TH} where TH is a vector of threshold 
%       percentages. This method compares the estimated onsets using
%       thresholding with each of the threshold values in TH and 
%       statistically determines the "best" onset from the list. If TH is
%       not specified, the following default vector of thresholds is used:
%       [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25].
%
%   `   7. {'multithresh2',TH} where TH is a vector of threshold 
%       percentages. This method compares the estimated onsets using
%       thresholding with each of the threshold values in TH before
%       estimating the "best" onset. If TH is not specified, the following 
%       default vector of thresholds is used: 5:5:60.
%
%       8. {'maxzeros'} estimates onset as the sample value which results 
%       in the IR having a maximum number of zeros inside the unit circle
%       when the IR is shifted to begin at that sample. The returned sample 
%       onset, O, is accurate to 1 sample. This option can result in a long
%       computation time.
%
%   [O,S] = ESTIMATEIRONSET(X,METHOD) when METHOD{1} is either 'grpdelay'
%   or 'phase' additionally returns the group delay or phase spectra,
%   respectively, where S has the same dimensions as X.
%
%   [O,N] = ESTIMATEIRONSET(X,{'maxzeros'}) additionally returns the
%   computed number of zeros inside the unit circle.
%
%   ___ = ESTIMATEIRONSET(...,'resample',F) resamples the IRs in X using a
%   resample factor of F. For example, if F is 2, then the IRs in X are
%   upsampled by a factor of 2. If F is 0.5, the IRs are downsampled by a
%   factor of 2. Returned onsets are in fractions of a sample at the
%   original sampling rate. Note that this option has an effect only when
%   METHOD{1} is one of the following: 'threshold', 'mpxc', 'mpxc2', 
%   'multithresh', 'multithresh2'.
%
%   Needs: Signal Processing Toolbox.
%
%   See also THRESHOLDIRS, GETGRPDELAY, MAKEMINPHASEIR.

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
%   Copyright (c) 2021 Princeton University
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

% Check input count
narginchk(1,4);

% Validate required input
validateattributes(inputIR,{'numeric'},{'2d','nonempty','nonnan',...
    'finite','real'},'estimateIROnset','X',1)
% If inputIR is a vector, force it to be a column vector.
inputIR = shiftdim(inputIR);

% Validate optional inputs
if nargin < 2
    METHOD = {'threshold',20};
else
    METHOD = varargin{1};
end
validateattributes(METHOD,{'cell'},{'nonempty'},'estimateIROnset',...
    'METHOD',2)
validateattributes(METHOD{1},{'char'},{'scalartext','nonempty'},...
    'estimateIROnset','the first element of METHOD',2)
numParams = length(METHOD);

indx = find(strcmpi(varargin,'resample'),1);
if indx
    resampFac = varargin{indx+1};
    if resampFac == 1
        resampIR = inputIR;
    else
        resampIR = resample(inputIR,resampFac,1);
    end
else
    resampFac = 1;
    resampIR = inputIR;
end

[irLen,numIRs] = size(inputIR); % This should be inputIR, not resampIR
switch lower(METHOD{1})
    case 'threshold'
        if numParams < 2
            thp = 20;
        else
            thp = METHOD{2};
            validateattributes(thp,{'numeric'},{'scalar','real',...
                'finite','nonnan','nonnegative','<=',100},...
                'estimateIROnset',['the THP value when METHOD Name is',...
                ' ''threshold'''],2)
        end
        
        onsetVal = thresholdIRs(resampIR,thp/100);
        onsetVal = onsetVal/resampFac;
    case 'grpdelay'
        if numParams > 1
            Fs = METHOD{2};
            validateattributes(Fs,{'numeric'},{'scalar','real','finite',...
                'nonnan','positive'},'estimateIROnset',['the value for',...
                ' Fs when METHOD Name is ''grpdelay'''],2)
        else
            error(['Sampling rate in Hz must be specified as an input',...
                ' when METHOD Name is ''grpdelay''.']);
        end
        
        if numParams < 3
            avgRange = [0,Fs/2];
        else
            avgRange = METHOD{3};
            validateattributes(avgRange,{'numeric'},{'vector','real',...
                'finite','nonnan','nonnegative','numel',2,'<=',Fs/2},...
                'estimateIROnset',['the values for [FL,FU] when METHOD',...
                ' Name is ''grpdelay'''],2)
        end
        
        [optOut,onsetVal] = getGrpDelay(inputIR,Fs,avgRange);
    case 'phase'
        if numParams > 1
            Fs = METHOD{2};
            validateattributes(Fs,{'numeric'},{'scalar','real','finite',...
                'nonnan','positive'},'estimateIROnset',['the value for',...
                ' Fs when METHOD Name is ''phase'''],2)
        else
            error(['Sampling rate in Hz must be specified as an input',...
                ' when METHOD Name is ''phase''.']);
        end
        
        fVec = getFreqVec(Fs,irLen);
        nyqFreq = Fs/2;
        if numParams < 3
            avgRange = [0,nyqFreq];
        else
            avgRange = METHOD{3};
            validateattributes(avgRange,{'numeric'},{'vector','real',...
                'finite','nonnan','nonnegative','numel',2,'<=',nyqFreq},...
                'estimateIROnset',['the values for [FL,FU] when METHOD',...
                ' Name is ''phase'''],2)
        end
        
        if avgRange(1) >= avgRange(2)
            error('FU must be greater than FL.')
        end
        
        [~,fLIndx] = min(abs(fVec-avgRange(1)));
        [~,fUIndx] = min(abs(fVec-avgRange(2)));
        if numParams > 3
            validateattributes(METHOD{4},{'char'},{'scalartext'},...
                'estimateIROnset',['the ''linearfit'' option when',...
                ' METHOD Name is ''phase'''],2)
            if strcmpi(METHOD{4},'linearfit')
                fitFlag = 1;
            else
                error(['Unrecognized parameter when METHOD Name is',...
                    ' ''phase''.']);
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
        
        nyqFreqIndx = ceil((irLen+1)/2);
        halfSpec = [zeros(1,numIRs);diag(1./(2*pi*fVec(2:nyqFreqIndx)))*...
            phaseSpec(2:nyqFreqIndx,:)];
        onsetVal = -Fs*mean(halfSpec(fLIndx:fUIndx,:),'omitnan');
        negHalfSpec = -flipud(halfSpec);
        if bitget(abs(irLen),1) % if irLen_orig is odd
            optOut = [halfSpec;negHalfSpec(1:(nyqFreqIndx-1),:)];
        else
            optOut = [halfSpec;negHalfSpec(2:(nyqFreqIndx-1),:)];
        end
    case 'mpxc'
        onsetVec = zeros(1,numIRs);
        minPhaseIR = makeMinPhaseIR(resampIR,'hilb');
        for ii = 1:numIRs
            [xc,lagVec] = xcoh(resampIR(:,ii),minPhaseIR(:,ii),...
                'timedomain');
            [~,lagIndex] = getInterpMax(abs(xc));
            x1 = floor(lagIndex);
            x2 = ceil(lagIndex);
            m = x2-lagIndex;
            onsetVec(ii) = (m*lagVec(x1))+((1-m)*lagVec(x2));
        end
        onsetVal = onsetVec;
        onsetVal = onsetVal/resampFac;
    case {'mpxc_robust','mpxc2'} % mpxc_robust for backwards compatibility
        switch numParams
            case 1
                error(['Sampling rate in Hz must be specified as an',...
                    ' input when METHOD Name is ''mpxc2''.']);
            case 2
                Fs = METHOD{2};
                FL = 500;
                FU = 1500;
                postMaxDur = 1;
                threshVal = 100;
                floorFlag = false;
            case 3
                Fs = METHOD{2};
                FL = min(METHOD{3});
                FU = max(METHOD{3});
                postMaxDur = 1;
                threshVal = 100;
                floorFlag = false;
            case 4
                Fs = METHOD{2};
                FL = min(METHOD{3});
                FU = max(METHOD{3});
                postMaxDur = METHOD{4};
                threshVal = 100;
                floorFlag = false;
            case 5
                Fs = METHOD{2};
                FL = min(METHOD{3});
                FU = max(METHOD{3});
                postMaxDur = METHOD{4};
                threshVal = METHOD{5};
                floorFlag = false;
            case 6
                Fs = METHOD{2};
                FL = min(METHOD{3});
                FU = max(METHOD{3});
                postMaxDur = METHOD{4};
                threshVal = METHOD{5};
                floorFlag = METHOD{6};
            otherwise
                error('Unrecognized inputs for METHOD ''mpxc2''.')
        end
        
        validateattributes(Fs,{'numeric'},{'scalar','real','finite',...
            'nonnan','positive'},'estimateIROnset',['the value for',...
            ' Fs when METHOD Name is ''mpxc2'''],2)
        validateattributes(FL,{'numeric'},{'scalar','real','finite',...
            'nonnan','nonnegative'},'estimateIROnset',['the value for',...
            ' FL when METHOD Name is ''mpxc2'''],2)
        validateattributes(FU,{'numeric'},{'scalar','real','finite',...
            'nonnan','positive'},'estimateIROnset',['the value for',...
            ' FU when METHOD Name is ''mpxc2'''],2)
        validateattributes(postMaxDur,{'numeric'},{'scalar','real',...
            'finite','nonnan','positive'},'estimateIROnset',['the',...
            ' value for postMaxDur when METHOD Name is ''mpxc2'''],2)
        validateattributes(threshVal,{'numeric'},{'scalar','real',...
            'finite','nonnan','nonnegative'},'estimateIROnset',['the',...
            ' value for threshVal when METHOD Name is ''mpxc2'''],2)
        validateattributes(floorFlag,{'logical'},{'scalar'},...
            'estimateIROnset',['the value for floorFlag when METHOD',...
            ' Name is ''mpxc2'''],2)
        
        Fs = resampFac*Fs;
        onsetVec = zeros(1,numIRs);
        postMaxLen = floor(postMaxDur*Fs/1000);
        inputIRLen = size(resampIR,1);
        maxIRLen = max([ceil(20*Fs/1000),inputIRLen]);
        for ii = 1:numIRs
            maxIndx = round(estimateIROnset(resampIR(:,ii),{'threshold',...
                threshVal}));
            winLenVec = maxIndx+postMaxLen;
            winIR = windowSignal(resampIR(:,ii),winLenVec,'wType',{'rc',...
                [0,postMaxLen/(2*winLenVec)]});
            winIR = padarray(winIR,[maxIRLen-winLenVec,0],0,'post');
%             minPhaseIR = makeMinPhaseIR(inputIR(:,ii),'hilb');
%             [xc,lagVec] = xcoh(inputIR(:,ii),minPhaseIR,'unweighted',...
%                 [FL,FU]/(Fs/2));
            minPhaseIR = makeMinPhaseIR(winIR,'hilb');
            [xc,lagVec] = xcoh(winIR,minPhaseIR,'unweighted',...
                [FL,FU]/(Fs/2));
%             [~,lagIndex] = max(abs(xc));
%             onsetVec(ii) = lagVec(lagIndex);
            [~,lagIndex] = getInterpMax(abs(xc));
            x1 = floor(lagIndex);
            x2 = ceil(lagIndex);
            m = x2-lagIndex;
            onsetVec(ii) = (m*lagVec(x1))+((1-m)*lagVec(x2));
        end
        
        onsetVec = onsetVec/resampFac;
        
        if floorFlag
            onsetVal = floor(onsetVec);
        else
            onsetVal = onsetVec;
        end
    case 'multithresh'
        if numParams < 2
            thpVec = [5;7.5;10;12.5;15;17.5;20;22.5;25];
        else
            thpVec = METHOD{2};
            validateattributes(thpVec,{'numeric'},{'vector','real',...
                'finite','nonnan','nonnegative','<=',100},...
                'estimateIROnset',['the TH value when METHOD Name is',...
                ' ''multithresh''.'],2)
        end
        
        % Compute onset using thresholding for each threshold in thpVec
        numTHPs = length(thpVec);
        delayMat = zeros(numTHPs,numIRs);
        for ii = 1:numTHPs
            delayMat(ii,:) = thresholdIRs(resampIR,thpVec(ii)/100);
        end
        
        % Check if all computed thresholds are equal
        delMatDiff = diff(delayMat);
        eqChkFlag = any(delMatDiff(:));
        if ~eqChkFlag
            onsetVal = delayMat(1,:);
        else
            finalDelayVec = zeros(1,numIRs);
            for ii = 1:numIRs
                currDelayVec = delayMat(:,ii);
                % Successively eliminate outlying delay values
                totDropped = 0;
                while (numTHPs-totDropped) >= 2
                    % Compute standard dev. excluding one delay at a time
                    numTHPsLeft = numTHPs-totDropped;
                    stDevVec = zeros(numTHPsLeft,1);
                    for jj = 1:numTHPsLeft
                        dropDelVec = currDelayVec;
                        dropDelVec(jj) = [];
                        stDevVec(jj) = std(dropDelVec);
                    end
                    
                    % If the standard devs. of the remaining delays are all 
                    % the same, exit the while loop
                    stdDiffVec = diff(stDevVec);
                    stdDiffChk = any(stdDiffVec);
                    if ~stdDiffChk
                        break
                    end
                    
                    % Find all indices of minimum standard deviation
                    drpIndxs = find(stDevVec == min(stDevVec));
                    drpCount = length(drpIndxs);
                    
                    % If the minimum standard deviation occurs at all 
                    % remaining positions, exit the while loop
                    if drpCount == numTHPsLeft
                        break
                    end
                    
                    % Drop elements
                    currDelayVec(drpIndxs) = [];
                    
                    totDropped = totDropped + drpCount;
                end
                finalDelayVec(ii) = currDelayVec(1);
            end
            onsetVal = finalDelayVec;
        end
        
        onsetVal = onsetVal/resampFac;
    case 'multithresh2'
        if numParams < 2
            thpVec = 5:5:60;
        else
            thpVec = METHOD{2};
            validateattributes(thpVec,{'numeric'},{'vector','real',...
                'finite','nonnan','nonnegative','<=',100},...
                'estimateIROnset',['the TH value when METHOD Name is',...
                ' ''multithresh2''.'],2)
        end
        
        % Compute onset using thresholding for each threshold in thpVec
        numTHPs = length(thpVec);
        delayMat = zeros(numTHPs,numIRs);
        for ii = 1:numTHPs
            delayMat(ii,:) = thresholdIRs(resampIR,thpVec(ii)/100);
        end
        
        delDiff = diff(delayMat);
        onsetVal = zeros(1,numIRs);
        for jj = 1:numIRs
            changeVec = [1;find(delDiff(:,jj));numTHPs];
            [~,maxChangeIndx] = max(diff(changeVec));
            onsetIndxL = changeVec(maxChangeIndx);
            onsetIndxU = changeVec(maxChangeIndx+1);
            thpValL = thpVec(onsetIndxL);
            thpValU = thpVec(onsetIndxU);
            thpValRange = thpValU-thpValL;
            tempOnset = delayMat(onsetIndxL,jj);
            if thpValL-thpValRange < 0
                onsetVal(1,jj) = tempOnset;
            else
                onsetVal(1,jj) = tempOnset-floor(thpValL/thpValRange);
            end
            onsetVal(1,jj) = clip(onsetVal(1,jj),1,tempOnset);
        end
        
        onsetVal = onsetVal/resampFac;
    case 'maxzeros'
        % First, get sample position of absolute maximum of each IR
        [~,absMaxVec] = max(abs(inputIR));
        
        onsetVec = absMaxVec;
        numZIn = zeros(1,numIRs);
        for ii = 1:numIRs
            cLB = 0;
            cUB = absMaxVec(ii);
            numShifts = cUB-cLB+1;
            stepSize = ceil(numShifts/10);
%             fprintf('Current bracket size: %d (LB: %d; UB: %d)\n',...
%                 numShifts,cLB,cUB);
            while numShifts > stepSize
                tSVec = round(linspace(cLB,cUB,stepSize));
                [tO,~] = getMZOnset(inputIR(1:absMaxVec(ii),ii),tSVec);
                tOIndx = find(tSVec == (tO-1),1,'last');
                if tOIndx == 1 || tOIndx == 6
                    cLB = tSVec(tOIndx);
                    cUB = tSVec(tOIndx);
                else
                    cLB = tSVec(tOIndx-1);
                    cUB = tSVec(tOIndx+1);
                end
                numShifts = cUB-cLB+1;
%                 fprintf('Current bracket size: %d (LB: %d; UB: %d)\n',...
%                     numShifts,cLB,cUB);
            end
            
            % Determine final shift vector
            shiftVec = cLB:cUB;
            [onsetVec(ii),numZIn(ii)] = getMZOnset(...
                inputIR(1:absMaxVec(ii),ii),shiftVec);
        end
        
        onsetVal = onsetVec;
        optOut = numZIn;
    otherwise
        error('Invalid METHOD Name specification.');
end

switch lower(METHOD{1})
    case {'grpdelay','phase','maxzeros'}
        nargoutchk(0,2);
        varargout{1} = onsetVal;
        varargout{2} = optOut;
    otherwise
        nargoutchk(0,1);
        varargout{1} = onsetVal;
end

end

function [O,maxN] = getMZOnset(X,S)
%GETMZONSET Estimate onset using 'max zeros' method.

% Compute number of zeros inside unit circle after each shift of the IR
numShifts = length(S);
numZIn = zeros(numShifts,1);
indx = 1;
for ii = S
    rX = circshift(X,-ii);
    rootVec = roots(rX);
    numZIn(indx) = numel(find(abs(rootVec) < 1));
    indx = indx + 1;
end

% Find position of max in numZIn
[maxN,maxNIndx] = max(numZIn);
O = S(maxNIndx)+1;

end
