function trueOnset = estimateTrueIROnset(h,approxOnset)
%ESTIMATETRUEIRONSET Estimate the true onset of an impulse response given 
%its approximate onset.
%   trueOnset = ESTIMATETRUEIRONSET(h,onset) estimates the true onset of 
%   IR in h given its approximate onset. If h is a matrix, the IRs must be 
%   stored as columns. The dimensions of approxOnset must be 1-by-n, where 
%   n is the number of columns in h. The dimensions of trueOnset are also 
%   1-by-n. The algorithm used is as follows and assumes that the specified 
%   approximate onset lies somewhere on the main impulse:
%       1. Locate first peak/notch at or after specified approximate onset.
%       2. Locate part of impulse up to first peak/notch where absolute
%       value of slope is maximum. It is assumed that this will occur
%       somewhere on the main impulse.
%       3. Estimate slope at location identified in step 2.
%       4. Find intercept (on time/sample axis) of the straight line with  
%       computed slope and passing through approximate onset point.
%       5. Intercept on time/sample axis is estimate of true onset.
%
%   Note: Unexpected results will be returned if the value of h at the 
%   specified approximate onset is not somewhere on the main impulse.
%
%   To compute approximate onsets, see ESTIMATEIRONSET.

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
%   Copyright (c) 2020 Princeton University
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

narginchk(2,2);

% Initialize required variables
[~,numIRs] = size(h);
trueOnset = zeros(1,numIRs);
h = [zeros(1,numIRs);h];

threshVec = (0.3:-0.05:0.05).';
% threshVec = 0;
numThreshs = length(threshVec);
for ii = 1:numIRs
    % Locate first peak/notch at or after specified approximate onset.
    maxAmp = max(abs(h(approxOnset(ii):end,ii)));
    [pMaxAmp,pIndxVec] = findpeaks(h(approxOnset(ii):end,ii));
    [nMaxAmp,nIndxVec] = findpeaks(-h(approxOnset(ii):end,ii));
    count = 1;
    loopFlag = true;
    while loopFlag && (count <= numThreshs)
        pIndx = find(pMaxAmp > threshVec(count)*maxAmp,1);
        nIndx = find(nMaxAmp > threshVec(count)*maxAmp,1);
        if ~isempty(pIndx) && ~isempty(nIndx)
            loopFlag = false;
            pnIndx = min([pIndxVec(pIndx),nIndxVec(nIndx)]);
        elseif isempty(nIndx) && ~isempty(pIndx)
            loopFlag = false;
            pnIndx = pIndxVec(pIndx);
        elseif isempty(pIndx) && ~isempty(nIndx)
            loopFlag = false;
            pnIndx = nIndxVec(nIndx);
        else
            count = count + 1;
        end
    end
    
    if loopFlag
        error('Unable to detect any peaks or notches.')
    end
    
    indxU = approxOnset(ii)+pnIndx-2; % Peak/notch location
    
    % Check indxU
    if indxU < 1
        error('Peak/notch detection algorithm failed.')
    end
    
    % Remove added zero
    hs = h(2:end,ii);
    
    % Estimate true onset
    switch indxU
        case {1,2}
            trueOnset(ii) = indxU;
        otherwise
            % Find peak or notch prior to the one located at indxU
            [~,priorPIndxVec] = findpeaks(hs(1:indxU));
            [~,priorNIndxVec] = findpeaks(-hs(1:indxU));
            indxL = max([max(priorPIndxVec),max(priorNIndxVec)]);
            if isempty(indxL)
                indxL = 1;
            end
            
            % Locate part of identified section of IR where absolute value
            % of slope is maximum, and find slope.
            slopeVec = diff(hs(indxL:indxU));
            absSlopeVec = abs(slopeVec);
            [~,maxSlopePos] = max(absSlopeVec);
            maxSlopeVal = slopeVec(maxSlopePos);
            
            % Estimate true onset
            indx = indxL+maxSlopePos;
            trueOnset(ii) = indx-(hs(indx)/maxSlopeVal);
    end
end

end
