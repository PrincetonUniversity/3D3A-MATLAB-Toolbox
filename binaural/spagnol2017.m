function DVF = spagnol2017(a,sPos,fVec,Fs,varargin)
%SPAGNOL2017 Distance variation function (DVF) by Spagnol et al. [1].
%   DVF = SPAGNOL2017(A,S,F,Fs) computes the DVF proposed by Spagnol et al.
%   [1] for a rigid sphere of radius A (specified in meters), sound source 
%   positions, S, in SOFA cartesian coordinates, vector of frequencies, F, 
%   in Hz, and sampling frequency, Fs, in Hz. S must be specified as an 
%   N-by-3 matrix. It is assumed that the ears are antipodal.
%
%   DVF = SPAGNOL2017(...,[EL;ER]) optionally specifies the positions of
%   the left and right "ears" on the surface of the sphere. EL and ER 
%   correspond to the left and right ear positions, respectively, each 
%   specified as the coordinate pair (az,el) where az and el are the 
%   azimuth and elevation, respectively, of the ear specified in degrees in 
%   SOFA spherical coordinates. 
%
%   DVF = SPAGNOL2017(...,[EL;ER],R) optionally specifies R as the source 
%   distance in meters. If R is not specified, R is computed directly from 
%   S. If R is specified, its value supersedes the source distance 
%   corresponding to the source positions given in S. R must be specified 
%   in meters. For a source infinitely far away, specify R as inf.

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

narginchk(4,6);

% Check inputs
validateattributes(a,{'numeric'},{'vector','nonempty','nonnan','finite',...
    'positive','real'},'spagnol2017','A',1);
validateattributes(sPos,{'numeric'},{'2d','nonempty','nonnan','finite',...
    'size',[NaN,3],'real'},'spagnol2017','S',2);
validateattributes(fVec,{'numeric'},{'vector','nonempty','nonnan',...
    'finite','nonnegative','real'},'spagnol2017','F',3);
validateattributes(Fs,{'numeric'},{'scalar','nonnan','finite',...
    'positive','real'},'spagnol2017','Fs',4);

% Get required variables and pre-process input data
a = shiftdim(a); % If a is a row vector, force it to be a column.
sDirs = sofaC2sofaS(sPos);
numDirs = size(sDirs,1); % numDirs = N in documentation.

if nargin < 6
    R = sDirs(:,3).';
else
    R = varargin{2};
    validateattributes(R,{'numeric'},{'vector','nonempty','nonnan'},...
        'spagnol2017','R',6);
    R = shiftdim(R).';
    
    if isscalar(R)
        R = R*ones(1,numDirs);
    else
        if size(R,1) ~= numDirs
            error('Invalid R specification.')
        end
    end
end
rho = (1./a)*R; % M-by-N matrix, where M is the number of 'a' values

if nargin < 5
    EL = [90,0];
    ER = [270,0];
else
    ePos = varargin{1};
    validateattributes(ePos,{'numeric'},{'2d','size',[2,2],'real'},...
        'spagnol2017','[EL;ER]',5);
    
    EL = ePos(1,:);
    ER = ePos(2,:);
end

% Compute angles of incidence
[xEL,yEL,zEL] = sofaS2sofaC(EL(1),EL(2),1);
[xER,yER,zER] = sofaS2sofaC(ER(1),ER(2),1);
S = sDirs(:,1:2);
[xS,yS,zS] = sofaS2sofaC(S(:,1),S(:,2),ones(numDirs,1));
thetaL = getCentralAngle([xS,yS,zS],repmat([xEL,yEL,zEL],numDirs,1));
thetaR = getCentralAngle([xS,yS,zS],repmat([xER,yER,zER],numDirs,1));
% Make compatible with dimensions of rho
thetaL = repmat(thetaL.',numAs,1);
thetaR = repmat(thetaR.',numAs,1);

% Main calculation

% Get vector of alpha (equivalent to theta) values for which coefficients 
% exist
aVec = getAlphaVec();

% Initialize variables
HL_DC = zeros(numAs,numDirs);
HR_DC = zeros(numAs,numDirs);
HL_HF = zeros(numAs,numDirs);
HR_HF = zeros(numAs,numDirs);
HL_Fc = zeros(numAs,numDirs);
HR_Fc = zeros(numAs,numDirs);
for ii = 1:numAs
    for jj = 1:numDirs
        % Perform calculations for left "ear".
        indx = find(aVec == thetaL(ii,jj),1);
        if ~isempty(indx)
            % Compute DC gain
            HL_DC(ii,jj) = computeDCGain(thetaL(ii,jj),rho(ii,jj));
            
            % Compute HF gain
            HL_HF(ii,jj) = computeHFGain(thetaL(ii,jj),rho(ii,jj));
            
            % Compute cut-off frequency
            HL_Fc(ii,jj) = computeCutoffFreq(thetaL(ii,jj),rho(ii,jj));
        else
            % Compute linearly-interpolated DC gain
            in1 = computeDCGain(10*floor(thetaL(ii,jj)/10),rho(ii,jj));
            in2 = computeDCGain(10*ceil(thetaL(ii,jj)/10),rho(ii,jj));
            HL_DC(ii,jj) = linearIterp(in1,in2,thetaL(ii,jj));
            
            % Compute linearly-interpolated HF gain
            in1 = computeHFGain(10*floor(thetaL(ii,jj)/10),rho(ii,jj));
            in2 = computeHFGain(10*ceil(thetaL(ii,jj)/10),rho(ii,jj));
            HL_HF(ii,jj) = linearIterp(in1,in2,thetaL(ii,jj));
            
            % Compute linearly-interpolated cut-off frequency
            in1 = computeCutoffFreq(10*floor(thetaL(ii,jj)/10),rho(ii,jj));
            in2 = computeCutoffFreq(10*ceil(thetaL(ii,jj)/10),rho(ii,jj));
            HL_Fc(ii,jj) = linearIterp(in1,in2,thetaL(ii,jj));
        end
        
        % Perform calculations for right "ear".
        indx = find(aVec == thetaR(ii,jj),1);
        if ~isempty(indx)
            % Compute DC gain
            HR_DC(ii,jj) = computeDCGain(thetaR(ii,jj),rho(ii,jj));
            
            % Compute HF gain
            HR_HF(ii,jj) = computeHFGain(thetaR(ii,jj),rho(ii,jj));
            
            % Compute cut-off frequency
            HR_Fc(ii,jj) = computeCutoffFreq(thetaR(ii,jj),rho(ii,jj));
        else
            % Compute linearly-interpolated DC gain
            in1 = computeDCGain(10*floor(thetaR(ii,jj)/10),rho(ii,jj));
            in2 = computeDCGain(10*ceil(thetaR(ii,jj)/10),rho(ii,jj));
            HR_DC(ii,jj) = linearIterp(in1,in2,thetaR(ii,jj));
            
            % Compute linearly-interpolated HF gain
            in1 = computeHFGain(10*floor(thetaR(ii,jj)/10),rho(ii,jj));
            in2 = computeHFGain(10*ceil(thetaR(ii,jj)/10),rho(ii,jj));
            HR_HF(ii,jj) = linearIterp(in1,in2,thetaR(ii,jj));
            
            % Compute linearly-interpolated cut-off frequency
            in1 = computeCutoffFreq(10*floor(thetaR(ii,jj)/10),rho(ii,jj));
            in2 = computeCutoffFreq(10*ceil(thetaR(ii,jj)/10),rho(ii,jj));
            HR_Fc(ii,jj) = linearIterp(in1,in2,thetaR(ii,jj));
        end
    end
end

% To be completed...

end

function aV = getAlphaVec()
%GETALPHAVEC Vector of alphas in Table 1 of Spagnol et al. [1].

narginchk(0,0);

aV = 0:10:180;

end

function dcVal = computeDCGain(alpha,rho)
%COMPUTEDCGAIN Compute DC gain.
%   D = COMPUTEDCGAIN(A,R) computes the DC gain, D, for a given angle of
%   incidence, A, and non-dimensional source distance, R, using Eq. (8) in
%   the work by Spagnol et al. [1].

narginchk(2,2);

dcVal = (p11(alpha)*rho+p21(alpha))/(rho^2+q11(alpha)*rho+q21(alpha));

end

function hfVal = computeHFGain(alpha,rho)
%COMPUTEHFGAIN Compute high-frequency gain.
%   HF = COMPUTEHFGAIN(A,R) computes the high-frequency gain, HF, for a 
%   given angle of incidence, A, and non-dimensional source distance, R, 
%   using Eq. (13) in the work by Spagnol et al. [1].

narginchk(2,2);

hfVal = (p12(alpha)*rho+p22(alpha))/(rho^2+q12(alpha)*rho+q22(alpha));

end

function fC = computeCutoffFreq(alpha,rho)
%COMPUTECUTOFFFREQ Compute cut-off frequency.
%   FC = COMPUTEHCOMPUTECUTOFFFREQFGAIN(A,R) computes the cut-off 
%   frequency, FC, in Hz for a given angle of incidence, A, and non-
%   dimensional source distance, R, using Eq. (14) in the work by Spagnol 
%   et al. [1].

narginchk(2,2);

fC = (p13(alpha)*rho^2+p23(alpha)*rho+p33(alpha))/(rho^2+q13(alpha)*rho+...
    q23(alpha));

end

function out = linearIterp(in1,in2,alpha)
%LINEARINTERP Perform linear interpolation.
%   O = LINEARINTERP(I1,I2,A) linearly interpolates between I1 and I2 based
%   on A and returns the interpolated value in O. See Eq. (15) in the work
%   by Spagnol et al. [1], for example.

narginchk(3,3);

out = ((ceil(alpha/10)-(alpha/10))*in1)+(((alpha/10)-floor(alpha/10))*in2);

end

function val = p11(alpha)
% P11 Value of p11 coefficient of RFA.
%   V = P11(A) returns the value, V, of the p11 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
p11Vec = [12.97,13.19,12.13,11.19,9.91,8.328,6.493,4.455,2.274,0.018,...
    -2.24,-4.43,-6.49,-8.34,-9.93,-11.3,-12.2,-12.8,-13];

indx = find(aVec == alpha,1);
val = p11Vec(indx);

end

function val = p21(alpha)
% P21 Value of p21 coefficient of RFA.
%   V = P21(A) returns the value, V, of the p21 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
p21Vec = [-9.69,234.2,-11.2,-9.03,-7.87,-7.42,-7.31,-7.28,-7.29,-7.48,...
    -8.04,-9.23,-11.6,-17.4,-48.4,9.149,1.905,-0.75,-1.32];

indx = find(aVec == alpha,1);
val = p21Vec(indx);

end

function val = q11(alpha)
% Q11 Value of q11 coefficient of RFA.
%   V = Q11(A) returns the value, V, of the q11 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
q11Vec = [-1.14,18.48,-1.25,-1.02,-0.83,-0.67,-0.5,-0.32,-0.11,-0.13,...
    0.395,0.699,1.084,1.757,4.764,-0.64,0.109,0.386,0.45];

indx = find(aVec == alpha,1);
val = q11Vec(indx);

end

function val = q21(alpha)
% Q21 Value of q21 coefficient of RFA.
%   V = Q21(A) returns the value, V, of the q21 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
q21Vec = [0.219,-8.5,0.346,0.336,0.379,0.421,0.423,0.382,0.314,0.24,...
    0.177,0.132,0.113,0.142,0.462,-0.14,-0.08,-0.06,-0.05];

indx = find(aVec == alpha,1);
val = q21Vec(indx);

end

function val = p12(alpha)
% P12 Value of p12 coefficient of RFA.
%   V = P12(A) returns the value, V, of the p12 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
p12Vec = [-4.39,-4.31,-4.18,-4.01,-3.87,-4.1,-3.87,-5.02,-6.72,-8.69,...
    -11.2,-12.1,-11.1,-11.1,-9.72,-8.42,-7.44,-6.78,-6.58];

indx = find(aVec == alpha,1);
val = p12Vec(indx);

end

function val = p22(alpha)
% P22 Value of p22 coefficient of RFA.
%   V = P22(A) returns the value, V, of the p22 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
p22Vec = [2.123,-2.78,4.224,3.039,-0.57,-34.7,3.271,0.023,-8.96,-58.4,...
    11.47,8.716,21.8,1.91,-0.04,-0.66,0.395,2.662,3.387];

indx = find(aVec == alpha,1);
val = p22Vec(indx);

end

function val = q12(alpha)
% Q12 Value of q12 coefficient of RFA.
%   V = Q12(A) returns the value, V, of the q12 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
q12Vec = [-0.55,0.59,-1.01,-0.56,0.665,11.39,-1.57,-0.87,0.37,5.446,...
    -1.13,-0.63,-2.01,0.15,0.243,0.147,-0.18,-0.67,-0.84];

indx = find(aVec == alpha,1);
val = q12Vec(indx);

end

function val = q22(alpha)
% Q22 Value of q22 coefficient of RFA.
%   V = Q22(A) returns the value, V, of the q22 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
q22Vec = [-0.06,-0.17,-0.02,-0.32,-1.13,-8.3,0.637,0.325,-0.08,-1.19,...
    0.103,-0.12,0.098,-0.4,-0.41,-0.34,-0.18,0.05,0.131];

indx = find(aVec == alpha,1);
val = q22Vec(indx);

end

function val = p13(alpha)
% P13 Value of p13 coefficient of RFA.
%   V = P13(A) returns the value, V, of the p13 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
p13Vec = [0.457,0.455,-0.87,0.465,0.494,0.549,0.663,0.691,3.507,-27.4,...
    6.371,7.032,7.092,7.463,7.453,8.101,8.702,8.925,9.317];

indx = find(aVec == alpha,1);
val = p13Vec(indx);

end

function val = p23(alpha)
% P23 Value of p23 coefficient of RFA.
%   V = P23(A) returns the value, V, of the p23 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
p23Vec = [-0.67,0.142,3404,-0.91,-0.67,-1.21,-1.76,4.655,55.09,10336,...
    1.735,40.88,23.86,102.8,-6.14,-18.1,-9.05,-9.03,-6.89];

indx = find(aVec == alpha,1);
val = p23Vec(indx);

end

function val = p33(alpha)
% P33 Value of p33 coefficient of RFA.
%   V = P33(A) returns the value, V, of the p33 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
p33Vec = [0.174,-0.11,-1699,0.437,0.658,2.02,6.815,0.614,589.3,16818,...
    -9.39,-44.1,-23.6,-92.3,-1.81,10.54,0.532,0.285,-2.08];

indx = find(aVec == alpha,1);
val = p33Vec(indx);

end

function val = q13(alpha)
% Q13 Value of q13 coefficient of RFA.
%   V = Q13(A) returns the value, V, of the q13 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
q13Vec = [-1.75,-0.01,7354,-2.18,-1.2,-1.59,-1.23,-0.89,29.23,1945,...
    -0.06,5.635,3.308,13.88,-0.88,-2.23,-0.96,-0.9,-0.57];

indx = find(aVec == alpha,1);
val = q13Vec(indx);

end

function val = q23(alpha)
% Q23 Value of q23 coefficient of RFA.
%   V = Q23(A) returns the value, V, of the q23 coefficient for a given
%   value of A based on the values provided by Spagnol et al. [1] in Table
%   1 of their paper where A corresponds to 'alpha' in the table.

narginchk(1,1);

aVec = getAlphaVec();
q23Vec = [0.699,-0.35,-5350,1.188,0.256,0.816,1.166,0.76,59.51,1707,...
    -1.12,-6.18,-3.39,-12.7,-0.19,1.295,-0.02,-0.08,-0.4];

indx = find(aVec == alpha,1);
val = q23Vec(indx);

end
