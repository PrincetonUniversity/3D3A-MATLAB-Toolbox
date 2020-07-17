function G = spagnol2017LFGain(a,sPos,varargin)
%SPAGNOL2017LFGAIN Approximate low-frequency gain for a rigid-sphere model.
%   G = SPAGNOL2017LFGAIN(A,S) computes an approx. low-frequency gain for a 
%   rigid sphere of radius A (specified in meters), and for sound source 
%   positions, S, specified in SOFA cartesian coordinates, using the
%   formula derived by Spagnol et al. [1]. S must be specified as an 
%   N-by-3 matrix. B will contain ILD values in dB. It is assumed that the
%   "ear" at which the gain is computed is at azimuth 90 and elevation 0.
%       If A is a scalar, B will be a row vector of length N.
%       If A is a vector of length M, B will be an M-by-N matrix.
%
%   G = SPAGNOL2017LFGAIN(...,E) optionally specifies the position of the
%   "ear" on the surface of the sphere at which the gain is computed. E
%   must be specified as a coordinate pair (az,el) where az and el 
%   correspond to the azimuth and elevation, respectively, of the ear 
%   specified in degrees in SOFA spherical coordinates. 
%
%   G = SPAGNOL2017LFGAIN(...,E,R) optionally specifies R as the source
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

narginchk(2,4);

% Check inputs
validateattributes(a,{'numeric'},{'vector','nonempty','nonnan','finite',...
    'positive','real'},'spagnol2017LFGain','A',1);
validateattributes(sPos,{'numeric'},{'2d','nonempty','nonnan','finite',...
    'size',[NaN,3],'real'},'spagnol2017LFGain','S',2);

% Get required variables and pre-process input data
a = shiftdim(a); % If a is a row vector, force it to be a column.
numAs = length(a);
sDirs = sofaC2sofaS(sPos);
numDirs = size(sDirs,1); % numDirs = N in documentation.

if nargin < 4
    R = sDirs(:,3).';
else
    R = varargin{2};
    validateattributes(R,{'numeric'},{'vector','nonempty','nonnan'},...
        'spagnol2017LFGain','R',4);
    R = shiftdim(R).';
    
    if isscalar(R)
        R = R*ones(1,numDirs);
    else
        if length(R) ~= numDirs
            error('Invalid R specification.')
        end
    end
end
rho = (1./a)*R; % M-by-N matrix, where M is the number of 'a' values

if nargin < 3
    E = [90,0];
else
    E = varargin{1};
    validateattributes(E,{'numeric'},{'vector','real'},...
        'spagnol2017LFGain','E',3);
end

% Compute angles of incidence
[xE,yE,zE] = sofaS2sofaC(E(1),E(2),1);
S = sDirs(:,1:2);
[xS,yS,zS] = sofaS2sofaC(S(:,1),S(:,2),ones(numDirs,1));
theta = getCentralAngle([xS,yS,zS],repmat([xE,yE,zE],numDirs,1));
% Make compatible with dimensions of rho
theta = repmat(theta.',numAs,1);

% Main calculation

% Get vector of alpha (equivalent to theta) values for which coefficients 
% exist
aVec = getAlphaVec();

% Initialize variables
G = zeros(numAs,numDirs);
for ii = 1:numAs
    for jj = 1:numDirs
        % Perform calculations for left "ear".
        indx = find(aVec == theta(ii,jj),1);
        if ~isempty(indx)
            % Compute DC gain
            G(ii,jj) = computeDCGain(theta(ii,jj),rho(ii,jj));
        else
            % Compute linearly-interpolated DC gain
            in1 = computeDCGain(10*floor(theta(ii,jj)/10),rho(ii,jj));
            in2 = computeDCGain(10*ceil(theta(ii,jj)/10),rho(ii,jj));
            G(ii,jj) = linearIterp(in1,in2,theta(ii,jj));
        end
    end
end

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
