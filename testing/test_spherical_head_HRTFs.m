% Testing spherical head HRTF computations

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

% Head parameters
a = 0.0875; % Head radius in m
eL = [90,0]; % Left ear direction
eR = [270,0]; % Right ear direction

% Source parameters
azVec = (0:50:350).';
elVec = 0;
rhoVec = [1.25;1.5;2;4;8;inf]; % See Fig. 3 in Duda and Martens [1].
rVec = a*rhoVec; % Source distance
rVecLen = length(rVec);

% DSP parameters
T = 0.005; % Duration of HRIRs in seconds
Fs = 44100;
method = {'exact',0};

%% Pre-compute source position and direction matrices

sourceDirs = makePosMat({azVec,elVec,1});
sourcePositions = sofaS2sofaC(sourceDirs);
sourcePositions = unique(round(sourcePositions,6),'rows','stable');
sourceDirs = sofaC2sofaS(sourcePositions);
% Update numDirs after removing duplicate directions
numDirs = size(sourceDirs,1);

%% Perform other pre-computation tasks

irLen = 2^nextpow2(T * Fs);
fVec = getFreqVec(Fs,irLen);
nyquistIndx = ceil((irLen+1)/2);

% Compute angles of incidence from source position data
fprintf('Computing angle(s) of incidence...')
S = sourceDirs(:,1:2);
numPos = size(S,1);
[xS,yS,zS] = sofaS2sofaC(S(:,1),S(:,2),ones(numPos,1));
% Compute left and right HRTFs if both ear positions are specified
[xEL,yEL,zEL] = sofaS2sofaC(eL(1),eL(2),1);
[xER,yER,zER] = sofaS2sofaC(eR(1),eR(2),1);
theta{1} = getCentralAngle([xS,yS,zS],repmat([xEL,yEL,zEL],numPos,1));
theta{2} = getCentralAngle([xS,yS,zS],repmat([xER,yER,zER],numPos,1));
fprintf('done!\n')

%% Compute spherical head HRTFs

numEars = numel(theta);
h = cell(numEars,1);
N = cell(numEars,1);
TH = cell(numEars,1);
hL = cell(rVecLen,1);
hR = cell(rVecLen,1);
NMatL = cell(rVecLen,1);
NMatR = cell(rVecLen,1);
pL1 = 100/rVecLen;
pL2 = pL1/numEars;
pL3 = pL2/numPos;
pL4 = pL3/nyquistIndx;
clearProgress = repmat('\b',1,6);
fprintf('Computing spherical head HRTFs...%5.1f%%',0);
for ll = 1:rVecLen
    for ii = 1:numEars
        H = zeros(nyquistIndx,numPos);
        N{ii,1} = zeros(nyquistIndx,numPos);
        TH{ii,1} = cell(nyquistIndx,numPos);
        for jj = 1:numPos
            for kk = 1:nyquistIndx
                % Progress
                fprintf(clearProgress);
                fprintf('%5.1f%%',((ll-1)*pL1 + (ii-1)*pL2 + ...
                    (jj-1)*pL3) + (kk-1)*pL4);
                [H(kk,jj),N{ii,1}(kk,jj),TH{ii,1}{kk,jj}] = ...
                    computeSphereHRTF(a,rVec(ll),theta{ii}(jj),fVec(kk),...
                    method);
            end
        end
        h{ii,1} = ifft(H,irLen,1,'symmetric');
    end
    hL{ll,1} = h{1,1};
    hR{ll,1} = h{2,1};
    NMatL{ll,1} = N{1,1};
    NMatR{ll,1} = N{2,1};
end
fprintf(clearProgress);
fprintf('%5.1f%%\n',100);

%% Initialize plot variables

irLen = size(hL{1,1},1);
fVec = getFreqVec(Fs,irLen);
nyqIndx = ceil((irLen+1)/2);
muVec = 2*pi*a*fVec/getSoundSpeed();
% Vector of angles of incidence w.r.t left ear
aoiVec = mod(sourceDirs(:,1)-eL(1),360); 
plotAngles = [0,90:10:180].';
numPlotAngles = length(plotAngles);

%% Generate plot similar to Fig. 2 in Duda and Martens [2].

plotIndxVec = zeros(numPlotAngles,1);
legendLabels = cell(numPlotAngles,1);
for ii = 1:numPlotAngles
    plotIndxVec(ii) = find(round(aoiVec) == plotAngles(ii));
    legendLabels{ii,1} = num2str(aoiVec(plotIndxVec(ii),1));
end
HL = getMagSpecdB(hL{rVecLen,1});

figure();
semilogx(muVec(1:(nyqIndx-1)),HL(1:(nyqIndx-1),plotIndxVec));
xlabel('Normalized frequency, $\mu$','Interpreter','latex')
ylabel('Response (dB)')
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

%% Generate plot similar to Fig. 3 in Duda and Martens [2].

HL = zeros(irLen,rVecLen*2);
legendLabels = cell(rVecLen*2,1);
for ii = 1:rVecLen
    HL(:,[ii,ii+rVecLen]) = getMagSpecdB(hL{ii,1}(:,[10,25]));
    legendLabels{ii,1} = ['\rho = ',num2str(rhoVec(ii)),...
        '; \theta = ',num2str(aoiVec(10))];
    legendLabels{ii+rVecLen,1} = ['\rho = ',num2str(rhoVec(ii)),...
        '; \theta = ',num2str(aoiVec(25))];
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
