% Testing spherical head HRTF computations
% Note: Cells have to be executed in sequence from top to bottom.
% Existing variables in the workspace will be cleared.

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

%   References: 
%       [1]. Duda and Martens (1998) - Range dependence of the response of 
%   a spherical head model.

% close all
clearvars
clc

% Head parameters
a = 0.09; % Head radius in m
eL = [100,0]; % Left ear direction
eR = [260,0]; % Right ear direction

% Source parameters
azVec = (0:10:350).';
elVec = 0;
rhoVec = [1.25;1.5;2;4;8;inf]; % See Fig. 3 in Duda and Martens [1].
rVec = a*rhoVec; % Source distance
rVecLen = length(rVec);

% DSP parameters
T = 0.01; % Duration of HRIRs in seconds
Fs = 48000;
causalFlag = true;
% Comment/uncomment the code below to test the different methods
% method = {'exact',inf};
% method = {'sridharchoueiri2019',3};
% method = {'cooperbauck1980',0.001};
method = {'dudamartens1998',0.001};
% method = {'fixedn',60};

%% Pre-compute source position and direction matrices

sourceDirs = makePosMat({azVec,elVec,1});
sourcePositions = sofaS2sofaC(sourceDirs);
[~,uniquePosIndxs,~] = unique(round(sourcePositions,6),'rows','stable');
sourcePositions = sourcePositions(uniquePosIndxs,:);
sourceDirs = sofaC2sofaS(sourcePositions);
% Update numDirs after removing duplicate directions
numDirs = size(sourceDirs,1);

%% Perform other pre-computation tasks

irLen = 2^nextpow2(T * Fs);
fVec = getFreqVec(Fs,irLen);
nyquistIndx = ceil((irLen+1)/2);

% Compute angles of incidence from source position data
fprintf('Computing angle(s) of incidence...')
[xS,yS,zS] = sofaS2sofaC(sourceDirs(:,1),sourceDirs(:,2),ones(numDirs,1));
% Compute left and right HRTFs if both ear positions are specified
[xEL,yEL,zEL] = sofaS2sofaC(eL(1),eL(2),1);
[xER,yER,zER] = sofaS2sofaC(eR(1),eR(2),1);
theta{1} = getCentralAngle([xS,yS,zS],repmat([xEL,yEL,zEL],numDirs,1));
theta{2} = getCentralAngle([xS,yS,zS],repmat([xER,yER,zER],numDirs,1));
fprintf('done!\n')

clear xS yS zS xEL yEL zEL xER yER zER

%% Compute spherical head HRTFs

numEars = numel(theta);
H = cell(numEars,rVecLen);
h = cell(numEars,1);
N = cell(numEars,1);
TH = cell(numEars,1);
hL = cell(rVecLen,1);
hR = cell(rVecLen,1);
NMatL = cell(rVecLen,1);
NMatR = cell(rVecLen,1);
threshL = cell(rVecLen,1);
threshR = cell(rVecLen,1);
pL1 = 100/rVecLen;
pL2 = pL1/numEars;
pL3 = pL2/numDirs;
pL4 = pL3/nyquistIndx;
clearProgress = repmat('\b',1,6);
fprintf('Computing spherical head HRTFs...%5.1f%%',0);
for ll = 1:rVecLen
    for ii = 1:numEars
        H{ii,ll} = zeros(nyquistIndx,numDirs);
        N{ii,1} = zeros(nyquistIndx,numDirs);
        TH{ii,1} = cell(nyquistIndx,numDirs);
        for jj = 1:numDirs
            for kk = 1:(nyquistIndx-1)
                % Progress
                fprintf(clearProgress);
                fprintf('%5.1f%%',((ll-1)*pL1 + (ii-1)*pL2 + ...
                    (jj-1)*pL3) + (kk-1)*pL4);
                [H{ii,ll}(kk,jj),N{ii,1}(kk,jj),TH{ii,1}{kk,jj}] = ...
                    computeSphereHRTF(a,rVec(ll),theta{ii}(jj),fVec(kk),...
                    method);
            end
            kk = nyquistIndx;
            % Progress
            fprintf(clearProgress);
            fprintf('%5.1f%%',((ll-1)*pL1 + (ii-1)*pL2 + ...
                (jj-1)*pL3) + (kk-1)*pL4);
            [H{ii,ll}(kk,jj),N{ii,1}(kk,jj),TH{ii,1}{kk,jj}] = ...
                computeSphereHRTF(a,rVec(ll),theta{ii}(jj),fVec(kk),...
                method);
            if mod(irLen,2) == 0 % irLen is even
                H{ii,ll}(kk,jj) = abs(H{ii,ll}(kk,jj));
            end
        end
        h{ii,1} = ifft(H{ii,ll},irLen,1,'symmetric');
    end
    hL{ll,1} = h{1,1};
    hR{ll,1} = h{2,1};
    if causalFlag && rVec(ll) == inf
        hL{ll,1} = shiftSignal(hL{ll,1},irLen/4);
        hR{ll,1} = shiftSignal(hR{ll,1},irLen/4);
    end
    NMatL{ll,1} = N{1,1};
    NMatR{ll,1} = N{2,1};
    threshL{ll,1} = TH{1,1};
    threshR{ll,1} = TH{2,1};
end
fprintf(clearProgress);
fprintf('%5.1f%%\n',100);

clear pL1 pL2 pL3 pL4 clearProgress ll ii jj kk h N TH

%% Initialize plot variables

irLen = size(hL{1,1},1);
fVec = getFreqVec(Fs,irLen);
nyqIndx = ceil((irLen+1)/2);
muVec = 2*pi*a*fVec/getSoundSpeed();
% Vector of angles of incidence w.r.t left ear
aoiVec = mod(round(sourceDirs(:,1)-eL(1),3),360); 
plotAngles = [0,90:10:180].';
numPlotAngles = length(plotAngles);

%% Generate plot similar to Fig. 2 in Duda and Martens [1].

% Get existing variables defined in workspace
varsbefore = who;

plotIndxVec = zeros(numPlotAngles,1);
legendLabels = cell(numPlotAngles,1);
for ii = 1:numPlotAngles
    plotIndxVec(ii) = find(round(aoiVec) == plotAngles(ii));
    legendLabels{ii,1} = num2str(aoiVec(plotIndxVec(ii),1));
end
HL = getMagSpecdB(hL{rVecLen,1});

figure();
semilogx(muVec(1:(nyqIndx)),HL(1:(nyqIndx),plotIndxVec));
xlabel('Normalized frequency, $\mu$','Interpreter','latex')
ylabel('Response (dB)','Interpreter','latex')
legend(legendLabels,'Location','southwest','Interpreter','latex')

ax = gca;
ax.XLim = [0.1,100];
ax.YLim = [-25,10];
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = ':';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'off';
ax.MinorGridLineStyle = ':';

% Clear variables generated only by this cell
varsafter = who;
newvars = setdiff(varsafter,varsbefore);
varstokeep = {}; % Variables specified here will not be cleared
varstoremove = setdiff(newvars,varstokeep);
clear(varstoremove{:});
clearvars varstokeep varstoremove varsafter newvars

%% Generate plot similar to Fig. 3 in Duda and Martens [1].

% Get existing variables defined in workspace
varsbefore = who;

HL = zeros(irLen,rVecLen*2);
legendLabels = cell(rVecLen*2,1);
for ii = 1:rVecLen
    HL(:,[ii,ii+rVecLen]) = getMagSpecdB(hL{ii,1}(:,[11,26]));
    legendLabels{ii,1} = ['\rho = ',num2str(rhoVec(ii)),...
        '; \theta = ',num2str(aoiVec(11))];
    legendLabels{ii+rVecLen,1} = ['\rho = ',num2str(rhoVec(ii)),...
        '; \theta = ',num2str(aoiVec(26))];
end

figure();
semilogx(muVec(1:(nyqIndx-1)),HL(1:(nyqIndx-1),:));
xlabel('Normalized frequency, $\mu$','Interpreter','latex')
ylabel('Response (dB)','Interpreter','latex')
legend(legendLabels,'Location','southwest')

ax = gca;
ax.XLim = [0.1,100];
ax.YLim = [-50,20];
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = ':';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'off';
ax.MinorGridLineStyle = ':';

% Clear variables generated only by this cell
varsafter = who;
newvars = setdiff(varsafter,varsbefore);
varstokeep = {}; % Variables specified here will not be cleared
varstoremove = setdiff(newvars,varstokeep);
clear(varstoremove{:});
clearvars varstokeep varstoremove varsafter newvars

%% Generate plot similar to Fig. 5 in Duda and Martens [1].
% Note: We use the same values of rho as before (and so some values don't
% match up with those in Fig. 5 in Duda and Martens [1]).

% Compute ILD
sPosIndx = 27; % This corresponds to an azimuth of 270 deg. (= 100 deg. in 
               % Duda and Martens [1]).
ILDCell = cell(rVecLen,1);
legendLabels = cell(rVecLen,1);
for ii = 1:rVecLen
    for jj = 1:nyquistIndx
        ILD = estimateILD(hL{ii,1}(:,sPosIndx),hR{ii,1}(:,sPosIndx));
        ILDCell{ii,1} = ILD.';
        legendLabels{ii,1} = ['\rho = ',num2str(rhoVec(ii))];
    end
end
ILDMat = (cell2mat(ILDCell)).';
% ILDMat_norm = ILDMat-repmat(ILDMat(1,:),irLen,1);
% ILDMat_norm2 = ILDMat_norm-repmat(ILDMat_norm(:,rVecLen),1,rVecLen);

clear ii jj sPosIndx ILDCell ILD

% Get existing variables defined in workspace
varsbefore = who;

figure();
semilogx(muVec(1:(nyqIndx-1)),ILDMat(1:(nyqIndx-1),:));
% semilogx(muVec(1:(nyqIndx-1)),ILDMat_norm2(1:(nyqIndx-1),:));
xlabel('Normalized frequency, $\mu$','Interpreter','latex')
ylabel('ILD (dB)','Interpreter','latex')
legend(legendLabels,'Location','southeast')

ax = gca;
ax.XLim = [0.1,100];
ax.YLim = [0,60];
% ax.YLim = [0,7];
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = ':';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'off';
ax.MinorGridLineStyle = ':';

% Clear variables generated only by this cell
varsafter = who;
newvars = setdiff(varsafter,varsbefore);
varstokeep = {}; % Variables specified here will not be cleared
varstoremove = setdiff(newvars,varstokeep);
clear(varstoremove{:});
clearvars varstokeep varstoremove varsafter newvars
clear legendLabels

%% Verify that returned order, N, is correct
%
% Method: Use returned order, N, and the 'fixedn' method to compute 
% thresholds which may then be compared to the input threshold. 
% Corresponding pairs of threshold values should match. This test is only 
% for 'cooperbauck1980' and 'dudamartens1998' methods.

thresh_test = cell(rVecLen,1);
pL1 = 100/rVecLen;
pL2 = pL1/numDirs;
pL3 = pL2/nyquistIndx;
clearProgress = repmat('\b',1,6);
fprintf('Performing verification...%5.1f%%',0);
for ll = 1:rVecLen
    thresh_test{ll,1} = cell(nyquistIndx,numDirs);
    for jj = 1:numDirs
        for kk = 1:nyquistIndx
            % Progress
            fprintf(clearProgress);
            fprintf('%5.1f%%',(ll-1)*pL1 + (jj-1)*pL2 + (kk-1)*pL3);
            [~,~,TH_test] = computeSphereHRTF(a,...
                rVec(ll),theta{1}(jj),fVec(kk),{'fixedn',...
                NMatL{ll,1}(kk,jj)+1}); % Plus 1 is needed and logical
            thresh_test{ll,1}{kk,jj} = threshL{ll,1}{kk,jj}-TH_test;
        end
    end
end
fprintf(clearProgress);
fprintf('%5.1f%%\n',100);

% Open thresh_test and see if the values are close to 0 (ignore any NaNs,
% if any, for mu = 0).

clear pL1 pL2 pL3 clearProgress ll jj kk TH_test
