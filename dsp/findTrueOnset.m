function trueOnset = findTrueOnset(h,fS,approxOnset)
%FINDTRUEONSET Compute the true onset of an impulse response given its
%approximate onset.
%   trueOnset = FINDTRUEONSET(h,fS,onset) computes the true onset of IR
%   in h given its approximate onset. If h is a matrix, the IRs must be 
%   stored as columns. The dimensions of approxOnset must be 1-by-n, where 
%   n is the number of columns in h. The dimensions of trueOnset are also 
%   1-by-n. fS is the sampling rate in Hz. The algorithm used is as
%   follows:
%       1. Locate first peak/notch after specified approximate onset.
%       2. Upsample input IRs and determine corresponding peak/notch sample 
%       location.
%       3. Find slope of attack at approx. onset sample.
%       4. Find intercept (on time/sample axis) of the straight line with  
%       computed slope and passing through approximate onset point.
%       5. Intercept on time/sample axis is the true onset.
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

narginchk(3,3);

% Upsample input IRs and locate first peak/notch after approximate onset.
[irLen,numIRs] = size(h);
u = 10; % Upsample factor.
uh = resample(h,u,1);
t = getTimeVec(fS,irLen);
ut = getTimeVec(u*fS,u*irLen);
firstPeakorNotch = zeros(1,numIRs);
uOnset = zeros(1,numIRs);
upsampPeakorNotch = zeros(1,numIRs);
meanSlope = zeros(1,numIRs);
trueOnset = zeros(1,numIRs);

signVec = sign(diff([zeros(1,numIRs);h])); % Sign of the derivative.
% Pre-padding with one zero accounts for the case when the first sample is
% the onset sample.
for ii = 1:numIRs
    jj = approxOnset(ii) + 1; % Plus 1 accounts for extra zero added above.
    
    % Compute first peak/notch after input onset sample.
    loopFlag = 1;
    while loopFlag && (jj <= (irLen+1))        
        if signVec(jj,ii) ~= signVec(jj-1,ii)
            firstPeakorNotch(ii) = jj - 1;
            loopFlag = 0;
        elseif signVec(jj,ii) ~= signVec(jj+1,ii)
            firstPeakorNotch(ii) = jj;
            loopFlag = 0;
        else
            jj = jj + 1;
        end
    end
    
    if loopFlag
        error(['Unable to detect first peak/notch after input sample ',...
            'for IR number %d in h',ii]);
    end
    
    % Find corresponding location in upsampled IR.
    [~,uOnset(ii)] = min(abs(ut-t(approxOnset(ii))));
    [~,upsampPeakorNotch(ii)] = min(abs(ut-t(firstPeakorNotch(ii))));
    
    % Find slope of attack
    lB = uOnset(ii)-round(u/2);
    % Calculation below accounts for case when onset = first peak/notch.
    uB = min([uOnset(ii)+round(u/2);upsampPeakorNotch(ii)]);
    attackSlice = abs(uh(lB:uB,ii));
    meanSlope(ii) = u*mean(gradient(attackSlice));
    
    % Compute true onset
    trueOnset(ii) = floor(approxOnset(ii)-(abs(h(approxOnset(ii),ii))/...
        meanSlope(ii)));
end

end

