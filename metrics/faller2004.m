function [tau_final,DL_final,fC,c12] = faller2004(bIn,fS,C0)
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
%   ___ = FALLER2004(A,FS,C0) optionally specifies the interaural coherence 
%   threshold (a value between 0 and 1) that is used to select "free-field"
%   ITD and ILD cues from the matrix of values that are computed. The
%   larger the C0, the closer the selected cues are to the free-field cues
%   (i.e. cues likely caused by reflections are ignored). However,
%   according to Faller and Merimaa [1], "there is a strong motivation to 
%   choose C0 as small as possible while still getting accurate enough 
%   ITD and/or ILD cues, because this will lead to the cues being selected 
%   more often, and consequently to a larger proportion of the ear input 
%   signals contributing to the localization." If C0 is not specified, a
%   default value of 0.99 is used.
%
%   Needs: LTFAT toolbox.

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

narginchk(2,3);

% Check inputs
validateattributes(bIn,{'double'},{'2d','nonempty','nonnan','finite',...
    'real','size',[NaN,2]},'faller2004','A',1);
validateattributes(fS,{'double'},{'scalar','nonempty','nonnan','finite',...
    'real','positive'},'faller2004','FS',2);

if nargin < 3
    C0 = 0.999;
end

validateattributes(C0,{'double'},{'scalar','nonempty','nonnan','finite',...
    'real','nonnegative','<=',1},'faller2004','C0',3);

% Auditory periphery modeling

% 1. Filter input binaural signal by Patterson's auditory filters with a
% uniform center frequency spacing of 1 ERB.
fC = getERBFreqVec(500,2000);
% fC = [500;2000]; % TODO: Make this an input
bIn_cbf = applyCBFilter(bIn,fS,fC); % bIn_cbf is a 2-element cell array.

% TODO: 1b. Add internal noise to model (see last para on p. 3077 of [1]).  

% 2. Apply envelope compression, half-wave rectification, squaring, and
% low-pass filtering to each signal. The resulting "x" signals correspond
% to "nerve firing densities". Apply these steps only to those signals
% corresponding to critical band center frequencies > 500.
x = bIn_cbf; % Initialize output
[lpf_num,lpf_den] = butter(4,425/(fS/2));
if any(fC > 500)
    for ii = 1:2
        cs = compressEnvelope(bIn_cbf{ii,1}(:,fC > 500),0.23);
        cs_hwr = rectifyInput(cs);
        cs_hwr_sq = cs_hwr.^2;
        x{ii,1}(:,fC > 500) = filter(lpf_num,lpf_den,cs_hwr_sq);
    end
end

% Binaural processor

% 1. Compute running cross-correlation, gamma
Tc = 0.01; % Time-constant
[sigLen,numCFs] = size(x{1,1});
gamma = cell(numCFs,1);
L1 = cell(numCFs,1); % Auto-correlation matrix for x{1,1}
L2 = cell(numCFs,1); % Auto-correlation matrix for x{2,1}
for ii = 1:numCFs
    [gamma{ii,1},L1{ii,1},L2{ii,1},lagVec] = runningXCorr(x{1,1}(:,ii),...
        x{2,1}(:,ii),fS,Tc);
end

% 2. Compute ITD matrix, tau, in seconds, and interaural coherence, c12
tau = zeros(sigLen,numCFs);
c12 = zeros(sigLen,numCFs);
for ii = 1:numCFs
    [c12(:,ii),lagIndx] = max(gamma{ii,1},[],2);
    tau(:,ii) = lagVec(lagIndx);
end
tau_s = tau/fS;

% 3. Compute ILD matrix, DL, in dB and total signal power contributing to
% ITD, ILD, IC cues
DL = zeros(sigLen,numCFs);
p = zeros(sigLen,numCFs);
for ii = 1:numCFs
    for jj = 1:sigLen
        DL(jj,ii) = mag2db(sqrt(L2{ii,1}(jj,tau(jj,ii)+max(lagVec)+1)./...
            L1{ii,1}(jj,tau(jj,ii)+max(lagVec)+1)));
        p(jj,ii) = L1{ii,1}(jj,tau(jj,ii)+max(lagVec)+1)+L2{ii,1}(jj,...
            tau(jj,ii)+max(lagVec)+1);       
    end
end

% Binaural cue extraction

% TODO: Allow for different C0 to be specified in each critical band

tau_final = cell(numCFs,1);
DL_final = cell(numCFs,1);
for ii = 1:numCFs
    validc12 = c12(:,ii) >= C0;
    tau_final{ii,1} = tau_s(validc12,ii);
    DL_final{ii,1} = DL(validc12,ii);
end

end
