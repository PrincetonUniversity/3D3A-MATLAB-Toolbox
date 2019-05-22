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
azVec = (0:10:350).';
elVec = 0;
sourceDirs = makePosMat({azVec,elVec,1});
rhoVec = [1.25;1.5;2;4;8;inf]; % See Fig. 3 in Duda and Martens [1].
rVec = a*rhoVec; % Source distance
rVecLen = length(rVec);

% DSP parameters
T = 0.02; % Duration of HRIRs in seconds
Fs = 96000;
N = 150; % Max. truncation order

%% Compute HRIRs

hL = cell(rVecLen,1);
hR = cell(rVecLen,1);
NMatL = cell(rVecLen,1);
NMatR = cell(rVecLen,1);
for ii = 1:rVecLen
    [hL{ii,1},hR{ii,1},NMatL{ii,1},NMatR{ii,1}] = ...
        computeSphericalHeadHRIRs(a,sourceDirs(:,1:2),'T',T,'fS',Fs,...
        'eL',eL,'eR',eR,'R',rVec(ii),'N',N,'autoselectN',true,'P',inf);
end

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
    plotIndxVec(ii) = find(aoiVec == plotAngles(ii));
    legendLabels{ii,1} = num2str(aoiVec(plotIndxVec(ii),1));
end
HL = getMagSpecdB(hL{rVecLen,1});

figure();
semilogx(muVec(1:(nyqIndx-1)),HL(1:(nyqIndx-1),plotIndxVec));
xlabel('Normalized frequency, $\mu$')
ylabel('Response (dB)')
legend(legendLabels,'Location','southwest')

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
xlabel('Normalized frequency, $\mu$')
ylabel('Response (dB)')
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
