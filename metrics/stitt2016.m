function [rE, rV] = stitt2016(sourcePositions, listenerPosition, sourceGains, alpha, timeDelays, sourceType, broadbandGains)
%STITT2016 Stitt's precedence-effect-based energy vector.
%   RE = STITT2016(RS,RL,G,ALPHA,T,TYPE,B) returns the K-by-3 matrix RE of
%   precedence-effect-based energy vectors (specified in Cartesian
%   coordinates) for K frequencies, as computed by Stitt et al. (2016).
%   
%   The inputs are as follows:
%    1) RS is an N-by-3 matrix of source positions (given in Cartesian
%       coordinates) for the N discrete sources;
%    2) RL is a 1-by-3 vector specifying the position of the listener;
%    3) G is a K-by-N matrix of complex-valued source gains, where each row
%       corresponds to a different frequency;
%    4) ALPHA is a scalar (0 <= ALPHA <= 1) that specifies the relative
%       importance of stationary signal energy to transient energy (see
%       paper for details);
%    5) T (optional) is an N-by-1 vector of additional time delays to be
%       added to each source signal;
%    6) TYPE (optional) is a string specifying the type of sources used
%       (either 'point' or 'plane');
%    7) B (optional) is a 1-by-N matrix of broadband source gains.
%
%   Note: if B is specified, then the precedence-effect weights are
%   computed using B as the source gains such that the resulting weights
%   are frequency-independent. Otherwise, the precedence-effect weights are
%   computed on a per-frequency basis.
%
%   [RE,RV] = STITT2016(...) additionally returns the K-by-3 matrix RV of
%   corresponding velocity vectors.
%
%   Note: if ALPHA is an array with NUMEL(ALPHA) > 1, then RE and RV are
%   returned as cell arrays with SIZE(RE) = SIZE(RV) = SIZE(ALPHA).
%
%   See also GERZON1992, TYLKA2017.

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
%   Copyright (c) 2019 Princeton University
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
%     [1] Stitt et al. (2016) Extended Energy Vector Prediction of
%         Ambisonically Reproduced Image Direction at Off-Center Listening
%         Positions.
%     [2] https://circlesounds.wordpress.com/matlab-code/

narginchk(4,7);

numSources = size(sourceGains,2);

if isempty(listenerPosition)
    listenerPosition = zeros(1,size(sourcePositions,2));
end

% Assumes point sources by default
if nargin < 6 || isempty(sourceType)
    sourceType = 'point';
end

if nargin < 7 || isempty(broadbandGains)
    broadbandFlag = false;
else
    broadbandFlag = true;
end

switch lower(sourceType)
    case {'point','point source','point-source'}
        % Compute source positions relative to the listener
        relativeSourcePositions = sourcePositions - ones(numSources,1)*listenerPosition;

        % Find distances of each source to the listener
        sourceDistances = sqrt(sum(relativeSourcePositions.^2,2));
        
        % Compute unit vectors pointing towards each source
        sourceDirections = normalizeVector(relativeSourcePositions,2);
        
        % Compute time-of-arrival for each source's signal
        arrivalTimes = sourceDistances/getSoundSpeed();
    case {'plane','plane wave','plane-wave'}
        % Set all distances to unity to ignore them in the calculation
        sourceDistances = ones(numSources,1);
        
        % Compute unit vectors pointing towards each source
        sourceDirections = normalizeVector(sourcePositions,2);
        
        % Compute time-of-arrival for each source's signal
        if norm(listenerPosition) == 0
            arrivalTimes = zeros(numSources,1);
        else
            additionalPathLengths = -dot(sourceDirections,ones(numSources,1)*listenerPosition,2);
            arrivalTimes = (additionalPathLengths - min(additionalPathLengths))/getSoundSpeed();
        end
    otherwise
        error('Unknown source type.');
end

if nargin < 5 || isempty(timeDelays)
    timeDelays = arrivalTimes;
else
    % Add in any user-specified additional delays
    timeDelays = timeDelays + arrivalTimes;
end

% sort signals
[sortedDelays, sortingOrder] = sort(timeDelays - min(timeDelays));
sortedDirections = sourceDirections(sortingOrder,:);
sortedDistances = sourceDistances(sortingOrder);
sortedGains = sourceGains(:,sortingOrder)*diag(1./sortedDistances);
if broadbandFlag
    sortedBroadbandGains = broadbandGains(:,sortingOrder)*diag(1./sortedDistances);
end

% remove zero-amplitude signals
if broadbandFlag
    nonzeroMask = sortedBroadbandGains~=0;
    sortedBroadbandGains = sortedBroadbandGains(nonzeroMask);
else
    nonzeroMask = any(sortedGains,1);
end
sortedGains = sortedGains(:,nonzeroMask);
sortedDirections = sortedDirections(nonzeroMask,:);
sortedDelays = sortedDelays(nonzeroMask);
sortedDelays = sortedDelays - min(sortedDelays);

% compute perceptual weights
if broadbandFlag
    % frequency-independent weights
    wTRThMat = ones(size(sortedGains,1),1) * computePrecedenceWeights(sortedBroadbandGains,sortedDirections,sortedDelays);
else
    % frequency-dependent weights
    wTRThMat = computePrecedenceWeights(sortedGains,sortedDirections,sortedDelays);
%     wTRThMat = zeros(size(sortedGains));
%     for kk = 1:size(sortedGains,1) % loop over frequency
%         wTRThMat(kk,:) = computePrecedenceWeights(sortedGains(kk,:),sortedDirections,sortedDelays);
%     end
end

rECell = cell(size(alpha));
rVCell = cell(size(alpha));
for aa = 1:numel(alpha)
    % Eq. (2)
    wMat = (1-alpha(aa))*wTRThMat + alpha(aa);
    
    % Eq. (1)
    sortedPowers = abs(wMat.*sortedGains).^2;
    rECell{aa} = diag(1./sum(sortedPowers,2))*sortedPowers*sortedDirections;
    
    % Apply to Gerzon's velocity vector
    weightedGains = wMat.*sortedGains;
    rVCell{aa} = real(diag(1./sum(weightedGains,2))*weightedGains*sortedDirections);
end

if numel(alpha) == 1
    rE = rECell{1};
    rV = rVCell{1};
else
    rE = rECell;
    rV = rVCell;
end

end

function wTRThMat = computePrecedenceWeights(sortedGains, sortedDirections, sortedDelays)
%computePrecedenceWeights Precedence-effect weights of Stitt et al. (2016).
%   W = computePrecedenceWeights(G, R, T) returns the K-by-N matrix W of
%   precedence-effect weights, for N discrete sources and K frequencies, as
%   computed by Stitt et al. (citation below).
%   
%   The three inputs are the K-by-N matrix G of source gains, N-by-3 matrix
%   R of source directions (specified as unit vectors in Cartesian
%   coordinates), and the length N vector T of time delays. Each input
%   array should be sorted such that n = 1:N lists the sources in order of
%   earliest arrival.
%   
%   W is K-by-N
%   G is K-by-N
%   R is N-by-3
%   T is N-by-1
%   
%   P. Stitt, S. Bertet, and M. van Walstijn, "Extended Energy Vector
%       Prediction of Ambisonically Reproduced Image Direction at Off-
%       Center Listening Positions," J. Audio Eng. Soc., vol. 64, no. 5,
%       pp. 299--310 (2016 May).

[kLen, numSources] = size(sortedGains);
sourceList = 1:numSources;
sigma = sqrt(sin(pi/45));

wTMat = zeros(numSources);
RMat = zeros(kLen,numSources);
wRMat = zeros(kLen,numSources);
wTh3DMat = zeros(numSources,numSources,kLen);
dThetaMat = zeros(numSources);
weightedGains = zeros(kLen,numSources);
RbarMat = zeros(kLen,numSources);
zetaMat = zeros(kLen,numSources);
wTRTh3DMat = zeros(numSources,numSources,kLen);
wTRThMat = zeros(kLen,numSources);
wTRThMat(:,sortedDelays == 0) = 1;
for ii = sourceList(sortedDelays > 0) % ii is the lagging source
    taui = 1000*sortedDelays(ii);
    
    % Eq. (6)
    Omegai = sourceList(sortedDelays < sortedDelays(ii));
%     Omegai = sourceList(sortedDelays <= sortedDelays(ii));
%     Omegai(Omegai==ii) = [];
    % Note: Omegai is the set of all leading sources
    
    tauj = 1000*sortedDelays(Omegai);
    % Eq. (5)
    wTMat(ii,Omegai) = erfc(1.09*(taui-tauj));
    
    % Eq. (8)
    RMat(:,ii) = 10*log10((abs(sortedGains(:,ii)).^2)./sum(abs(sortedGains(:,Omegai)).^2,2));
    
    % Eq. (7)
    temp = 10.^(RMat(:,ii)/20) - 4.5;
    wRMat(:,ii) = 0.5 + sign(temp).*gammainc((abs(temp)/3.5).^10,1/10,'lower')/2;
    
    % Eq. (10)
    dThetaMat(ii,Omegai) = (sortedDirections(Omegai,2) - sortedDirections(ii,2)).';
    % Note: assumes y axis (second coordinate) is the interaural axis
    
    % Eq. (9)
    wTh3DMat(ii,Omegai,:) = (wTRThMat(:,Omegai)*diag(exp(-0.5*(dThetaMat(ii,Omegai)/sigma).^2))).';
    
    % Eq. (11)
    wTRTemp = wRMat(:,ii)*(1 - wTMat(ii,Omegai)) + ones(kLen,1)*wTMat(ii,Omegai);
    
    % Eq. (12)
    wTRTh3DMat(ii,Omegai,:) = (1 - shiftdim(wTh3DMat(ii,Omegai,:),1)).*wTRTemp.' + shiftdim(wTh3DMat(ii,Omegai,:),1);
    
    % Eq. (16)
    weightedGains(:,Omegai) = wTRThMat(:,Omegai).*sortedGains(:,Omegai);
    RbarMat(:,Omegai) = 10*log10(diag(1./sum(abs(weightedGains(:,Omegai)).^2,2))*(abs(weightedGains(:,Omegai)).^2));
    
    % Eq. (15)
    zetaMat(:,Omegai) = 10.^(RbarMat(:,Omegai)/20);
    
    % Eq. (14)
    wTRThMat(:,ii) = sum(zetaMat(:,Omegai).*(shiftdim(wTRTh3DMat(ii,Omegai,:),1).'),2)./sum(zetaMat(:,Omegai),2);
end

wTRThMat(isnan(wTRThMat)) = 0;

end