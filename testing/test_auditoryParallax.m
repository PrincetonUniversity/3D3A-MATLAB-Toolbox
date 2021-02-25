%% Testing auditory parallax computations
%
% *NOTE*: 
%   1. Cells have to be executed in sequence from top to bottom.
%   2. Existing variables in the workspace will be cleared, open figures 
%   will be closed, and the Command Window will be cleared. To avoid this,
%   do NOT execute the code in this cell, but the code may not run
%   correctly.
%
% Tested functions:
%   1. auditoryParallax

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

close all
clearvars
clc

%% Define input positions, PI

az = 0:5:355;
el = 0;
r = 0.18;
sDirs = makePosMat({az;el;r});
pIn = sofaS2sofaC(sDirs);

clearvars -except pIn sDirs r

%% Define ear position, E, and compute auditory parallax

eDir = [90,0,0.09];
ePos = sofaS2sofaC(eDir);
R = 2*r; % When equal to 'r', pIH_dirs and sDirs should match.

[pO_inf,pO_R] = auditoryParallax(pIn,ePos,R);

%% Visualize results

pO_inf_dirs = sofaC2sofaS(pO_inf);
pIH_dirs = sofaC2sofaS(pO_R);

figure();
plot(sDirs(:,1),[sDirs(:,1),pO_inf_dirs(:,1),pIH_dirs(:,1)]);
