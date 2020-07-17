function G = computeSphereLFGain(a,sPos,varargin)
%COMPUTESPHERELFGAIN Exact low-frequency gain for a rigid-sphere model.
%   G = COMPUTESPHERELFGAIN(A,S) computes exact low-frequency gain for a 
%   rigid sphere of radius A (specified in meters), and for sound source 
%   positions, S, specified in SOFA cartesian coordinates, using the
%   formula derived by Sridhar and Choueiri. S must be specified as an 
%   N-by-3 matrix. B will contain ILD values in dB. It is assumed that the
%   "ear" at which the gain is computed is at azimuth 90 and elevation 0.
%       If A is a scalar, B will be a row vector of length N.
%       If A is a vector of length M, B will be an M-by-N matrix.
%
%   G = COMPUTESPHERELFGAIN(...,E) optionally specifies the position of the
%   "ear" on the surface of the sphere at which the gain is computed. E
%   must be specified as a coordinate pair (az,el) where az and el 
%   correspond to the azimuth and elevation, respectively, of the ear 
%   specified in degrees in SOFA spherical coordinates. 
%
%   G = COMPUTESPHERELFGAIN(...,E,R) optionally specifies R as the source
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
    'positive','real'},'computeSphereLFGain','A',1);
validateattributes(sPos,{'numeric'},{'2d','nonempty','nonnan','finite',...
    'size',[NaN,3],'real'},'computeSphereLFGain','S',2);

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
        'computeSphereLFGain','R',4);
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
        'computeSphereLFGain','E',3);
end

% Compute angles of incidence
[xE,yE,zE] = sofaS2sofaC(E(1),E(2),1);
S = sDirs(:,1:2);
[xS,yS,zS] = sofaS2sofaC(S(:,1),S(:,2),ones(numDirs,1));
theta = getCentralAngle([xS,yS,zS],repmat([xE,yE,zE],numDirs,1));
% Make compatible with dimensions of rho
theta = repmat(theta.',numAs,1);

% DC gain calculation
g = sqrt(rho.^2 - 2*rho.*cosd(theta) + 1);
G = (2*rho./g)-(rho.*log((g+1-(rho.*cosd(theta)))./(rho.*(1-...
    cosd(theta)))));
zeroIndx = find(theta == 0,1);
if ~isempty(zeroIndx)
    G(:,zeroIndx) = (2*rho(:,zeroIndx)./g(:,zeroIndx))-(rho(:,...
        zeroIndx).*log(rho(:,zeroIndx)./(rho(:,zeroIndx)-1)));
end
% if theta == 0
%     G = (2*rho./g)-(rho.*log(rho./(rho-1)));
% else
%     G = (2*rho./g)-(rho.*log((g+1-(rho.*cosd(theta)))./(rho.*(1-...
%         cosd(theta)))));
% end

end
