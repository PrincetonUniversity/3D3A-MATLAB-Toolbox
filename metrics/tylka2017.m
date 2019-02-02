function [rPAvg, rPVec, fc] = tylka2017(a, Fs, VECTOR, varargin)
%TYLKA2017 Tylka's precedence-effect localization model for ambisonics.
%   [RAVG,RVEC,FC] = TYLKA2017(A,FS) computes predicted localization
%   vectors RVEC at each frequency FC, as well as a stimulus-weighted
%   average RAVG taken over frequency, given input ambisonics signals A
%   specified at a sampling rate FS.
%
%   [...] = TYLKA2017(A,FS,TYPE) computes localization vectors of TYPE
%   'energy', 'velocity', or 'combined'.
%
%   [...] = TYLKA2017(A,FS,TYPE,Name1,Value1,Name2,Value2,...) specifies
%   optional comma-separated pairs of Name,Value arguments, where Name is
%   the argument name and Value is the corresponding value. Name must
%   appear inside single quotes (' '). You can specify several name and
%   value pair arguments in any order as Name1,Value1,...,NameN,ValueN.
%   Valid Name,Value arguments are as follows:
%
%   'Band Average'          Two values: bandwidth (either 'erb' or a scalar
%                           for a fractional octave bandwidth) and a
%                           frequency range pair.
%
%   'Stimulus'              Stimulus signal column vector.
%
%   'AmbNorm'               String specifying the ambisonics normalization
%                           of signals A.
%
%   'Crossover'             Crossover frequency for the 'combined' vector.
%
%   'Energy Parameters'     Three values: alpha parameter, plane-wave grid
%                           file, and wavelet parameters for the energy
%                           vector.
%
%   'Velocity Parameters'   Same as 'Energy Parameters' for the velocity
%                           vector.
%
%   See also STITT2016, COMPUTEBANDAVG.

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

%   References:
%     [1] Tylka and Choueiri (2017) Models for evaluating navigational
%         techniques for higher-order ambisonics.
%     [2] Stitt et al. (2016) Extended Energy Vector Prediction of
%         Ambisonically Reproduced Image Direction at Off-Center Listening
%         Positions.

if nargin < 3 || isempty(VECTOR)
    VECTOR = 'combined';
end

% Specify band averaging parameters
indx = find(strcmpi(varargin,'Band Average'));
if indx
    bandWidth = varargin{indx+1}; % fractional-octave or ERB
    avgRange = varargin{indx+2};
else
    bandWidth = 'erb'; % fractional-octave or ERB
    avgRange = [20 20000];
end

% Specifiy stimulus signal
indx = find(strcmpi(varargin,'Stimulus'));
if indx
    stim = varargin{indx+1};
    stimLen = size(stim,1);
    stimPS = abs(fft(stim,stimLen,1)).^2;
else
    stimLen = 2048; % arbitrarily chosen length
    stimPS = ones(stimLen,1);
end
[stimWeights, fc] = computeBandAvg(stimPS,getFreqVec(Fs,stimLen),bandWidth,avgRange,Fs);
stimWeights = stimWeights/sum(stimWeights); % normalize

% Specifiy ambisonics normalization
indx = find(strcmpi(varargin,'AmbNorm'));
if indx
    ambNorm = varargin{indx+1};
else
    ambNorm = 'N3D';
end

% Determine crossover frequency band
indx = find(strcmpi(varargin,'Crossover'));
if indx
    xoFreq = varargin{indx+1};
else
    xoFreq = 700;
end
if xoFreq < fc(1)
    xoIndx = 0;
elseif xoFreq > fc(end)
    xoIndx = length(fc);
else
    [~, xoIndx] = findNearest(fc,xoFreq);
end

% Specific energy vector parameters
indx = find(strcmpi(varargin,'Energy Parameters'));
if indx
    alpha_E = varargin{indx+1};
    pwGridFile_E = varargin{indx+2};
    prepParams_E = varargin{indx+3};
else
    alpha_E = 0.7;
    pwGridFile_E = 'fliege_36';
    prepParams_E = {1e-3, -30, 500, 0.1};
end

% Specific velocity vector parameters
indx = find(strcmpi(varargin,'Velocity Parameters'));
if indx
    alpha_V = varargin{indx+1};
    pwGridFile_V = varargin{indx+2};
    prepParams_V = varargin{indx+3};
else
    alpha_V = 0.95;
    pwGridFile_V = 'fliege_9';
    prepParams_V = {2e-3, -8, 500, 0.1};
end

switch lower(VECTOR)
    case 'energy'
        rPVec = compute_rP(a, Fs, VECTOR, alpha_E, pwGridFile_E, prepParams_E, bandWidth, avgRange, ambNorm);
    case 'velocity'
        rPVec = compute_rP(a, Fs, VECTOR, alpha_V, pwGridFile_V, prepParams_V, bandWidth, avgRange, ambNorm);
    case 'combined'
        rPE = compute_rP(a, Fs,  'energy' , alpha_E, pwGridFile_E, prepParams_E, bandWidth, avgRange, ambNorm);
        rPV = compute_rP(a, Fs, 'velocity', alpha_V, pwGridFile_V, prepParams_V, bandWidth, avgRange, ambNorm);
        rPVec = combineVectors(rPV,rPE,xoIndx);
end
rPAvg = stimWeights.'*rPVec;

end

function rP = compute_rP(a, Fs, VECTOR, alpha, pwGridFile, prepParams, bandWidth, avgRange, ambNorm)
% Compute single vector spectrum

[pwGrid, wQList] = loadGridFile(pwGridFile);
muQList = a2mu(a,pwGrid,ambNorm,false)*diag(wQList); % Compute plane-wave decomposition

% Break-up plane-wave IRs into wavelets
[fullpwGrid, pwSourceGains, timeDelays] = prepIRs(pwGrid, muQList, Fs, prepParams);

% Compute average amplitude in specified frequency bands
[pwSourceGains, ~] = computeBandAvg(abs(pwSourceGains),getFreqVec(Fs, size(a,1)),bandWidth,avgRange,Fs);

switch lower(VECTOR)
    case 'energy'
        [rP, ~] = stitt2016(fullpwGrid, [], pwSourceGains, alpha, timeDelays, 'plane');
    case 'velocity'
        [~, rP] = stitt2016(fullpwGrid, [], pwSourceGains, alpha, timeDelays, 'plane');
end

end

function rC = combineVectors(rV,rE,xoIndx)
% Combine vector spectra

switch xoIndx
    case 0
        rC = rE;
    case size(rV,1)
        rC = rV;
    otherwise
        % scale velocity vectors to match magnitude of energy vectors at crossover
        normFactor = norm(rE(xoIndx,:))/norm(rV(xoIndx,:));
        rC = [normFactor*rV(1:xoIndx,:); rE((xoIndx+1):end,:)];
end

end

function [allPositions, allGains, allDelays] = prepIRs(sourcePositions, sourceIRs, Fs, params)
% Break up source IRs into wavelets
%   [RSo, G, D] = prepIRs(RSi, IR)
%
%   The outputs are as follows:
%    1) RSo is an No-by-3 matrix of source positions (in Cartesian
%       coordinates) for the No discrete sources;
%    2) G is a K-by-No matrix of source gains;
%    3) D is an No-by-1 vector of additional time delays corresponding to
%       each source signal;
%
%   The inputs are as follows:
%    1) RSi is an Ni-by-3 matrix of source positions (in Cartesian
%       coordinates) for the Ni discrete sources;
%    2) IR is a K-by-Ni matrix of source signals;

if nargin < 4 || isempty(params)
    % Arbitrarily chosen default parameters
    winTailLen = round(Fs*1e-3); % 1 ms
    peakThresh = db2mag(-18); % -18 dB (12.6%)
    hpfCutoff = 500;
    onsetThresh = 0.1;
else
    winTailLen = round(Fs*params{1});
    peakThresh = db2mag(params{2});
    hpfCutoff = params{3};
    onsetThresh = params{4};
end

[allIRs, allSampDelays, indexVec] = isolateWavelets(sourceIRs, winTailLen, peakThresh, hpfCutoff/(Fs/2), onsetThresh);
allPositions = sourcePositions(indexVec,:);
allGains = getPotential(allIRs,[],1);
allDelays = allSampDelays/Fs;

end

function [Y, D, L] = isolateWavelets(X, winTailLen, peakThresh, hpfCutoff, onsetThresh)
% Split IRs into wavelets; find delays

if nargin < 5
    onsetThresh = 0.1;
end

if nargin < 4 || isempty(hpfCutoff)
    hpfFlag = false;
else
    hpfFlag = true;
    [hpfB,hpfA] = butter(4,hpfCutoff,'high');
end

[IRLen, Ni] = size(X);
peakPosCell = cell(Ni,1);

No = 0;
globalPeak = max(max(abs(X)));
for ii = 1:Ni
    if hpfFlag % High-pass filter each IR to find peaks
        Xii = filter(hpfB, hpfA, X(:,ii));
    else
        Xii = X(:,ii);
    end
    
    if max(abs(Xii)) >= peakThresh*globalPeak
        [~, peakPosCell{ii}] = findpeaks(abs(Xii)/globalPeak, ...
            'MinPeakHeight', peakThresh, 'MinPeakDistance', winTailLen);
    else
        peakPosCell{ii} = [];
    end
    
    if isempty(peakPosCell{ii})
        No = No + 1;
    else
        No = No + length(peakPosCell{ii});
    end
end

Y = zeros(IRLen,No);
L = zeros(No,1);

indx = 1;
for ii = 1:Ni
    if isempty(peakPosCell{ii}) % No prominent peak; copy entire signal
        L(indx) = ii;
        Y(:,indx) = X(:,ii);
        indx = indx + 1;
    else % At least one prominent peak; window and copy each portion
        numPeaks = length(peakPosCell{ii});
        for jj = 1:numPeaks
            L(indx) = ii;
            
            winStart = max(peakPosCell{ii}(jj) - winTailLen, 1);
            if jj == numPeaks
                winEnd = IRLen;
            else
                winEnd = peakPosCell{ii}(jj+1);
            end
            
            % Window should be at least two tail-lengths long
            winLen = winEnd - winStart + 1;
            if winLen < 2*winTailLen
                winEnd = min(winStart + 2*winTailLen - 1, IRLen);
                winLen = winEnd - winStart + 1;
            end
            taperFrac = min(2*winTailLen/winLen, 1);
            winRange = winStart:winEnd;
            
            Y(winRange,indx) = tukeywin(winLen,taperFrac).*X(winRange,ii);
            
            indx = indx + 1;
        end
    end
end

if hpfFlag % High-pass filter each IR before finding onsets
    D = (thresholdIRs(filter(hpfB,hpfA,Y,[],1), onsetThresh)-1).';
else
    D = (thresholdIRs(Y, onsetThresh)-1).';
end

end