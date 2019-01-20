function [tau_final,DL_final,fC,c12,p_final] = faller2004(bIn,fS,fC,C0,...
    Tc,compFac,NOISE,SPL)
%FALLER2004 Implementation of Faller and Merimaa's binaural localization
%model.
%   B = FALLER2004(A,FS) applies Faller and Merimaa's binaural localization
%   model to the input binaural signal A, which must be an N-by-2 matrix,
%   where N is the length of the binaural signal. The sampling rate, FS,
%   must also be specified in Hz. B is a 2-element cell array where the
%   elements correspond to critical band center frequencies of 500 and 2000
%   Hz. Each element contains a vector of ITD values in seconds and should
%   contain ITDs of all sources in the sound field. However, binning of ITD 
%   values corresponding to different sources needs to be performed
%   externally.
%
%   B = FALLER2004(...,FC) optionally specifies a vector of critical band
%   center frequencies at which to estimate the binaural cues. Specify FC
%   as [] for default of FC = [500,2000].
%
%   B = FALLER2004(...,C0) optionally specifies the interaural coherence 
%   threshold (a value between 0 and 1) that is used to select "free-field"
%   ITD and ILD cues from the matrix of values that are computed. The
%   larger the C0, the closer the selected cues are to the free-field cues
%   (i.e. cues likely caused by reflections are ignored). However,
%   according to Faller and Merimaa [1], "there is a strong motivation to 
%   choose C0 as small as possible while still getting accurate enough 
%   ITD and/or ILD cues, because this will lead to the cues being selected 
%   more often, and consequently to a larger proportion of the ear input 
%   signals contributing to the localization." If C0 is not specified, a
%   default value of 0.98 is used.
%       If C0 is a scalar, the same value is used at all critical band
%       center frequencies. 
%       If C0 is a vector, it must have length equal to 2 if FC is not
%       specified, or the same length as the specified FC.
%   Specify C0 as [] for default of C0 = 0.98.
%
%   B = FALLER2004(...,TC) optionally specifies the time constant, in
%   seconds, for the running cross-correlation performed in the binaural
%   processing stage. The default value is 0.01 (i.e. 10 ms). Specify TC as
%   [] for default.
%
%   B = FALLER2004(...,CF) optionally specifies the compression factor
%   for compressing the envelope. CF can take values between 0 and 1. It is
%   set to 0.23 by default.
%
%   B = FALLER2004(...,NOISE) optionally specifies whether or not to add
%   independent Gaussian noise to the model to describe the limited 
%   accuracy of the auditory system. The options are: 'nonoise' (default)
%   and 'addnoise'. Specify NOISE as [] for default.
%
%   B = FALLER2004(...,SPL) specifies the average SPL of the input signal.
%   This value is useful only when 'addnoise' is specified as the previous
%   input. If 'addnoise' is specified but an SPL value is not, then SPL is
%   assumed to be 90 dB.
%
%   [B,D] = FALLER2004(...) additionally returns ILD values in dB. The
%   storage format for D is identical to that of B. Note that, according to 
%   Faller and Merimaa [1], due to envelope compression, these ILD 
%   estimates will be smaller than the level differences between the ear 
%   input signals by a factor of 4.3478 (equal to 1/0.23, where 0.23 is the 
%   compression factor).
%
%   [B,D,FC] = FALLER2004(...) additionally returns the vector of center
%   frequencies in Hz.
%
%   [B,D,FC,C] = FALLER2004(...) additionally returns the interaural
%   coherence matrix C, which is an N-by-2 matrix with N being the number
%   of samples in the input binaural signal, and 2 being the number of
%   critical band center frequencies.
%
%   [B,D,FC,C,P] = FALLER2004(...) additionally returns the total signal
%   energy, P, of the signals that contribute to the estimated ITD, ILD and
%   IC cues at each time index and for each critical band center frequency.

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

%   Ref:
%       [1]. Faller and Merimaa (2004) - Source localization in complex 
%       listening situations/ Selection of binaural cues based on 
%       interaural coherence.

narginchk(2,8);

% Check inputs
validateattributes(bIn,{'double'},{'2d','nonempty','nonnan','finite',...
    'real','size',[NaN,2]},'faller2004','A',1);
bLen = size(bIn,1);

validateattributes(fS,{'double'},{'scalar','nonempty','nonnan','finite',...
    'real','positive'},'faller2004','FS',2);

if nargin < 8
    SPL = 90;
end

if nargin < 7 || isempty(NOISE)
    NOISE = 'nonoise';
end

if nargin < 6 || isempty(compFac)
    compFac = 0.23;
end

if nargin < 5 || isempty(Tc)
    Tc = 0.01;
end
validateattributes(Tc,{'double'},{'scalar','nonempty','nonnan','finite',...
    'real','positive'},'faller2004','TC',5);

if nargin < 4 || isempty(C0)
    C0 = 0.98;
end
validateattributes(C0,{'double'},{'vector','nonempty','nonnan','finite',...
    'real','nonnegative','<=',1},'faller2004','C0',4);

if nargin < 3 || isempty(fC)
    fC = [500; 2000];
end
validateattributes(fC,{'double'},{'vector','nonempty','nonnan','finite',...
    'real','nonnegative','<=',fS/2},'faller2004','FC',3);
fC = shiftdim(fC);
numCFs = length(fC);

if ~isscalar(C0)
    if length(C0) ~= numCFs
        error('Expected length of C0 is either 1 or %d samples.',...
            length(fC))
    else
        C0 = shiftdim(C0);
    end
else
    C0 = C0*ones(numCFs,1);
end

% Auditory periphery modeling

% 1. Filter input binaural signal by Patterson's auditory filters with
% center frequencies specified in fC.
bIn_cbf = applyCBFilter(bIn,fS,fC); % bIn_cbf is a 2-element cell array.

% 1b. Optionally add independent Gaussian noise
energyL = sum(abs(bIn_cbf{1,1}).^2);
energyR = sum(abs(bIn_cbf{2,1}).^2);
avgEnergy = 0.5*(energyL+energyR);
if strcmpi(NOISE,'addnoise')
    noiseL = randn(bLen,1);
    noiseR = randn(bLen,1);
    noise_cbf = applyCBFilter([noiseL,noiseR],fS,fC);
    
    % Hearing threshold according to ISO 389 1975
    % Borrowed from MATLAB code written by Faller and Merimaa
    hthr = [
    %    Hz    dB SPL (re: 2e-5 Pa)
        125     47.9
        250     28.3
        500     12.6
        1000    6.8
        1500    6.7
        2000    7.8
        3000    7.6
        4000    8.7
        6000    11.9
        8000    11.6
    ];
    thr = zeros(numCFs,1);
    for ii = 1:numCFs
        [~,indx] = min(abs(hthr(:,1)-fC(ii)));
        thr(ii) = hthr(indx,2);
        
        gain = sqrt(avgEnergy/sum(noise_cbf{1,1}(:,ii).^2)).*...
            10^((thr(ii)-SPL)/20);
        bIn_cbf{1,1}(:,ii) = bIn_cbf{1,1}(:,ii) + ...
            gain.*noise_cbf{1,1}(:,ii);
        
        gain = sqrt(avgEnergy/sum(noise_cbf{2,1}(:,ii).^2)).*...
            10^((thr(ii)-SPL)/20);
        bIn_cbf{2,1}(:,ii) = bIn_cbf{2,1}(:,ii) + ...
            gain.*noise_cbf{2,1}(:,ii);
    end
end

% 2. Apply envelope compression, half-wave rectification, squaring, and
% low-pass filtering to each signal. The resulting "x" signals correspond
% to "nerve firing densities". Apply these steps only to those signals
% corresponding to critical band center frequencies >= 500 Hz.
x = bIn_cbf; % Initialize output
[lpf_num,lpf_den] = butter(4,425/(fS/2));
f_cutoff = 500; % Cut-off frequency beyond which the following is applied.
if any(fC >= f_cutoff)
    for ii = 1:2 % Each channel
        cs = compressEnvelope(bIn_cbf{ii,1}(:,fC >= f_cutoff),compFac);
        cs_hwr = rectifyInput(cs);
        cs_hwr_sq = cs_hwr.^2;
        x{ii,1}(:,fC >= f_cutoff) = filter(lpf_num,lpf_den,cs_hwr_sq);
    end
end

% Binaural processor

% 1. Compute running cross-correlation, gamma
gamma = cell(numCFs,1);
L1 = cell(numCFs,1); % Auto-correlation matrix for x{1,1}
L2 = cell(numCFs,1); % Auto-correlation matrix for x{2,1}
for ii = 1:numCFs
    [gamma{ii,1},L1{ii,1},L2{ii,1},lagVec] = runningXCorr(x{1,1}(:,ii),...
        x{2,1}(:,ii),fS,Tc);
end

% 2. Compute ITD matrix, tau, in seconds, and interaural coherence, c12
tau = zeros(bLen,numCFs);
c12 = zeros(bLen,numCFs);
for ii = 1:numCFs
    [c12(:,ii),lagIndx] = max(gamma{ii,1},[],2);
    tau(:,ii) = lagVec(lagIndx);
end
tau_s = tau/fS; % Negative sign -> sources on the left have negative ITD

% 3. Compute ILD matrix, DL, in dB and total signal power contributing to
% ITD, ILD, IC cues
DL = zeros(bLen,numCFs);
p = zeros(bLen,numCFs);
for ii = 1:numCFs
    for jj = 1:bLen
        % Sources on the left should have negative DL
        DL(jj,ii) = mag2db(sqrt(L2{ii,1}(jj,tau(jj,ii)+max(lagVec)+1)./...
            L1{ii,1}(jj,tau(jj,ii)+max(lagVec)+1)));
        p(jj,ii) = L1{ii,1}(jj,tau(jj,ii)+max(lagVec)+1)+L2{ii,1}(jj,...
            tau(jj,ii)+max(lagVec)+1);       
    end
end

% Binaural cue extraction

tau_final = cell(numCFs,1);
DL_final = cell(numCFs,1);
p_final = cell(numCFs,1);

% Compute offset for returning data (mainly to avoid initial time segment
% where SNR may be very low)
% offsetL = find(cumsum(abs(bIn_cbf{1,1}).^2) > 0.0001*avgEnergy,1);
% offsetR = find(cumsum(abs(bIn_cbf{2,1}).^2) > 0.0001*avgEnergy,1);
% offset = max([offsetL,offsetR]);
offset = 1;
for ii = 1:numCFs
    indxList = double(c12(offset:end,ii) >= C0(ii));
    indxDiff = diff(indxList);
    lIndx = find(indxDiff,1);
    if isempty(lIndx)
        lIndx = 1;
        hIndx = bLen;
    else
        hIndx = find(indxDiff == -1,1);
    end
    tau_final{ii,1} = tau_s(lIndx:hIndx,ii);
    DL_final{ii,1} = DL(lIndx:hIndx,ii);
    p_final{ii,1} = p(lIndx:hIndx,ii);
end

end
