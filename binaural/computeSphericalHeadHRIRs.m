function [hL,hR,NMatL,NMatR,varargout] = computeSphericalHeadHRIRs(a,S,...
    varargin)
%COMPUTESPHERICALHEADHRIRS HRIRs of a rigid spherical head.
%   [hL,hR,NMatL,NMatR] = COMPUTESPHERICALHEADHRIRS(A,S) returns HRIRs for  
%   a spherical head of radius A and for a source at direction S. The HRIRs 
%   have a duration of 0.005 s at a sampling rate of 44100 Hz. The left and 
%   right "ear" positions are at [90,0,A] and [270,0,A], respectively, when 
%   specified in SOFA spherical coordinates. The source distance is set to 
%   infinity. The infinite sum in the analytical formulation (see refs. 
%   [1-3]) is truncated to 41 terms by default. The head radius should be 
%   specified in metres. The source direction should be specified in SOFA 
%   spherical coordinates as [az,el], in degrees. If S is a matrix, HRIRs
%   are computed for each source direction. S must have dimensions P-by-2,
%   where P is the number of directions. Also returned are NMatL and NMatR,
%   matrices of truncation orders used at each frequency up to the Nyquist
%   frequency for computing HRTFs corresponding to the left (NMatL) and 
%   right (NMatR) "ears". The values in NMatL and NMatR are useful when the 
%   optional flag 'autoselectN' is set to true (see below).
%
%   ___ = COMPUTESPHERICALHEADHRIRS(...,Name1,Value1,...) allows optional
%   specification of the following Name-Value pairs:
%       1. 'T',T - duration of the HRIRs in seconds.
%       2. 'fS',fS - sampling rate in Hz.
%       3. 'eL',eL - left ear position specified in SOFA
%       spherical coordinates as [az,el], in degrees.
%       4. 'eR',eR - right ear position specified in SOFA
%       spherical coordinates as [az,el], in degrees.
%       5. 'R',R - source distance in meters (must be > A).
%       6. 'N',N - truncation order for truncating the infinite sum to
%       have N+1 terms, or the maximum truncation order if 'autoselectN' is
%       set to true (see option 8 below).
%       7. 'makecausal',c - flag to specify whether or not to force output 
%       IRs to be causal by applying an artificial delay. c can be true or 
%       false (default).
%       8. 'autoselectN',f - flag to specify whether or not to autoselect
%       the truncation order for each frequency subject to a maximum value
%       of N, where N is input by the user separately (see option 6 above).
%       The autoselection algorithm chooses the highest order that produces
%       a change in the computed HRTF at that frequency. f can be true or 
%       false (default).
%       9. 'P',p - precision to use when autoselectN is set to true. p must
%       be a non-negative integer. The specified value of p is used to 
%       round the differences between successive partial sums when 
%       determing N. If p is specified as inf (default), no rounding
%       (beyond intrinsic numerical precision constraints) is performed.
%
%   [___,dTL,dTR] = COMPUTESPHERICALHEADHRIRS(...) optionally returns
%   matrices of threshold values (one for each ear) used in the algorithm 
%   by Duda and Martens [1]. The matrices have N+1 columns and give the
%   threshold value for each successive value of N as a function of 
%   frequency.
%
%   Needs: Symbolic Math Toolbox, MATLAB R2014b or later.

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
%       [1]. Duda and Martens (1998) - Range dependence of the response of 
%   a spherical head model.
%       [2]. Cooper (1982) - Calculator Program for Head-Related Transfer 
%   Function.
%       [3]. Morse and Ingard (1968) - Theoretical Acoustics.

narginchk(2,20);

% Parse and verify inputs
inputs = parseComputeSphericalHeadHRIRsInputs(a,S,varargin);

% Extract parsed inputs
a = inputs.a;
S = inputs.S;
T = inputs.T;
fS = inputs.fS;
eL = inputs.eL;
eR = inputs.eR;
R = inputs.R;
N = inputs.N;
cFlag = inputs.makecausal;
nFlag = inputs.autoselectN;
pres = inputs.P;

% Perform additional checks on inputs
if fS < 1000
    warning(['Input sampling rate, %d, is low. Note that sampling',...
        ' rate must be specified in Hz.'],fS)
end
if eL(1) < 0
    error('Invalid az specification for eL. az must be within [0,360).')
end
if eR(1) < 0
    error('Invalid az specification for eR. az must be within [0,360).')
end
if eL(2) > 90
    error('Invalid el specification for eL. el must be within [-90,90].')
end
if eR(2) > 90
    error('Invalid el specification for eR. el must be within [-90,90].')
end

% Initialize variables
irLen = round(T*fS);
fVec = getFreqVec(fS,irLen);
nyquistIndx = ceil((irLen+1)/2);
c = getSoundSpeed();
mu = fVec(1:nyquistIndx)*2*pi*a/c;
numPos = size(S,1);
[xS,yS,zS] = sofaS2sofaC(S(:,1),S(:,2),ones(numPos,1));
[xEL,yEL,zEL] = sofaS2sofaC(eL(1),eL(2),1);
[xER,yER,zER] = sofaS2sofaC(eR(1),eR(2),1);
mVec = (0:N).';
thetaL = zeros(numPos,1);
thetaR = zeros(numPos,1);
hrtfL = ones(nyquistIndx,numPos);
hrtfR = ones(nyquistIndx,numPos);
PL = zeros(length(mVec),numPos);
PR = zeros(length(mVec),numPos);
psiL = zeros(nyquistIndx,numPos);
psiR = zeros(nyquistIndx,numPos);
NMatL = zeros(nyquistIndx,numPos);
NMatR = zeros(nyquistIndx,numPos);

% Compute term in solution associated with source direction dependence
clearProgress = repmat('\b',1,6);
progressLimit1 = round(100*(1-(numPos/nyquistIndx)));
fprintf('Computing spherical head HRTFs...%5.1f%%',0);
for ii = 1:numPos
    % Progress
    fprintf(clearProgress);
    fprintf('%5.1f%%',(ii-1)*progressLimit1/numPos);
    thetaL(ii) = getCentralAngle([xS(ii),yS(ii),zS(ii)],[xEL,yEL,zEL]);
    thetaR(ii) = getCentralAngle([xS(ii),yS(ii),zS(ii)],[xER,yER,zER]);
    PL(:,ii) = legendreP(mVec,cosd(thetaL(ii)));
    PR(:,ii) = legendreP(mVec,cosd(thetaR(ii)));
end

dTL = cell(nyquistIndx,1);
dTR = cell(nyquistIndx,1);
if R == inf % Source at infinity. Following Cooper [2].
    for ii = 2:nyquistIndx
        % Progress
        fprintf(clearProgress);
        fprintf('%5.1f%%',progressLimit1+((ii-1)*(100-progressLimit1)/...
            nyquistIndx));
        
        dh = dSphericalHankelH(mVec,1,mu(ii));
        % conj so that negative phase = delay
        psi = diag(conj((2*mVec+1).*((-1i).^(mVec-1))./dh));
        psiPL = psi*PL;
        psiPR = psi*PR;
        
        dTL{ii,1} = (psiPL./cumsum(psiPL)).';
        dTR{ii,1} = (psiPR./cumsum(psiPR)).';
        if nFlag % Auto-select N if true
            if pres ~= inf
                changeVecL = round(diff(cumsum(psiPL)),pres);
                changeVecR = round(diff(cumsum(psiPR)),pres);
            else
                changeVecL = diff(cumsum(psiPL));
                changeVecR = diff(cumsum(psiPR));
            end
            for jj = 1:numPos
                nonNaNIndxs = find(~isnan(changeVecL(:,jj)));
                lastNonZeroChange = find(changeVecL(nonNaNIndxs,jj),1,...
                    'last');
                if isempty(lastNonZeroChange)
                    maxNL = 0;
                else
                    % Note: maxNL is the order, and thus is calculated as 
                    % the index of the last non-zero change and not the 
                    % index 1 after that, as might be expected.
                    maxNL = nonNaNIndxs(lastNonZeroChange);
                end
                NMatL(ii,jj) = min([maxNL,N]);
                
                nonNaNIndxs = find(~isnan(changeVecR(:,jj)));
                lastNonZeroChange = find(changeVecR(nonNaNIndxs,jj),1,...
                    'last');
                if isempty(lastNonZeroChange)
                    maxNR = 0;
                else
                    maxNR = nonNaNIndxs(lastNonZeroChange);
                end
                NMatR(ii,jj) = min([maxNR,N]);
            end
        else
            NMatL(ii,:) = N*ones(1,numPos);
            NMatR(ii,:) = N*ones(1,numPos);
        end
        for jj = 1:numPos
            psiL(ii,jj) = sum(psiPL(1:(NMatL(ii,jj)+1),jj),1);
            psiR(ii,jj) = sum(psiPR(1:(NMatR(ii,jj)+1),jj),1);
        end
        if cFlag
            % extra exp term to make IRs causal
            del = mu(ii)*(c/a)*(irLen/fS)*0.5;
            hrtfL(ii,:) = (1/(mu(ii)^2))*exp(-1i*del)*psiL(ii,:);
            hrtfR(ii,:) = (1/(mu(ii)^2))*exp(-1i*del)*psiR(ii,:);
        else
            hrtfL(ii,:) = (1/(mu(ii)^2))*psiL(ii,:);
            hrtfR(ii,:) = (1/(mu(ii)^2))*psiR(ii,:);
        end
    end
else
    rho = R/a;
    for ii = 2:nyquistIndx
        % Progress
        fprintf(clearProgress);
        fprintf('%5.1f%%',progressLimit1+((ii-1)*(100-progressLimit1)/...
            nyquistIndx));
        
        h = sphericalHankelH(mVec,1,mu(ii)*rho);
        dh = dSphericalHankelH(mVec,1,mu(ii));
        psi = diag(conj(((2*mVec)+1).*(h./dh)));
        psiPL = psi*PL;
        psiPR = psi*PR;
        
        dTL{ii,1} = (psiPL./cumsum(psiPL)).';
        dTR{ii,1} = (psiPR./cumsum(psiPR)).';
        if nFlag
            if pres ~= inf
                changeVecL = round(diff(cumsum(psiPL)),pres);
                changeVecR = round(diff(cumsum(psiPR)),pres);
            else
                changeVecL = diff(cumsum(psiPL));
                changeVecR = diff(cumsum(psiPR));
            end
            for jj = 1:numPos
                nonNaNIndxs = find(~isnan(changeVecL(:,jj)));
                lastNonZeroChange = find(changeVecL(nonNaNIndxs,jj),1,...
                    'last');
                if isempty(lastNonZeroChange)
                    maxNL = 0;
                else
                    maxNL = nonNaNIndxs(lastNonZeroChange);
                end
                NMatL(ii,jj) = min([maxNL,N]);
                
                nonNaNIndxs = find(~isnan(changeVecR(:,jj)));
                lastNonZeroChange = find(changeVecR(nonNaNIndxs,jj),1,...
                    'last');
                if isempty(lastNonZeroChange)
                    maxNR = 0;
                else
                    maxNR = nonNaNIndxs(lastNonZeroChange);
                end
                NMatR(ii,jj) = min([maxNR,N]);
            end
        else
            NMatL(ii,:) = N*ones(1,numPos);
            NMatR(ii,:) = N*ones(1,numPos);
        end
        for jj = 1:numPos
            psiL(ii,jj) = sum(psiPL(1:(NMatL(ii,jj)+1),jj),1);
            psiR(ii,jj) = sum(psiPR(1:(NMatR(ii,jj)+1),jj),1);
        end
        if cFlag
            hrtfL(ii,:) = (rho/mu(ii))*exp(-1i*mu(ii)*rho)*psiL(ii,:);
            hrtfR(ii,:) = (rho/mu(ii))*exp(-1i*mu(ii)*rho)*psiR(ii,:);
        else
            hrtfL(ii,:) = (rho/mu(ii))*psiL(ii,:);
            hrtfR(ii,:) = (rho/mu(ii))*psiR(ii,:);
        end
    end
end

hL = ifft(hrtfL,irLen,1,'symmetric');
hR = ifft(hrtfR,irLen,1,'symmetric');

if nargout > 4
    varargout{1} = dTL;
    varargout{2} = dTR;
end

fprintf(clearProgress);
fprintf('%5.1f%%\n',100);

end

function inputs = parseComputeSphericalHeadHRIRsInputs(a,S,opts)
%PARSECOMPUTESPHERICALHEADHRIRSINPUTS Parse and verify inputs for the 
%computeSphericalHeadHRIRs function.

p = inputParser;

% Required inputs
addRequired(p,'a',@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','nonnan','finite','real','positive'},...
    'computeSphericalHeadHRIRs','A',1));
addRequired(p,'S',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite','real','ncols',2,'>=',-90,'<',360},...
    'computeSphericalHeadHRIRs','S',2));

% Optional inputs
addParameter(p,'T',0.005,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','nonnan','finite','real','positive'},...
    'computeSphericalHeadHRIRs','T'));
addParameter(p,'fS',44100,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','nonnan','finite','real','positive'},...
    'computeSphericalHeadHRIRs','fS'));
addParameter(p,'eL',[90,0],@(x)validateattributes(x,{'double'},...
    {'vector','nonempty','nonnan','finite','real','numel',2,'>=',-90,...
    '<',360},'computeSphericalHeadHRIRs','eL'));
addParameter(p,'eR',[270,0],@(x)validateattributes(x,{'double'},...
    {'vector','nonempty','nonnan','finite','real','numel',2,'>=',-90,...
    '<',360},'computeSphericalHeadHRIRs','eR'));
addParameter(p,'R',inf,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','nonnan','real','positive','>',a},...
    'computeSphericalHeadHRIRs','R'));
addParameter(p,'N',40,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','nonnan','finite','nonnegative','integer'},...
    'computeSphericalHeadHRIRs','N'));
addParameter(p,'makecausal',false,@(x)validateattributes(x,{'logical'},...
    {'nonempty','nonnan'},'computeSphericalHeadHRIRs','makecausal, c'));
addParameter(p,'autoselectN',false,@(x)validateattributes(x,...
    {'logical'},{'nonempty','nonnan'},'computeSphericalHeadHRIRs',...
    'autoselectN, f'));
addParameter(p,'P',inf,@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','nonnegative'},'computeSphericalHeadHRIRs','P'));

p.CaseSensitive = false;
p.FunctionName = 'computeSphericalHeadHRIRs';

parse(p,a,S,opts{:});

inputs = p.Results;

end
