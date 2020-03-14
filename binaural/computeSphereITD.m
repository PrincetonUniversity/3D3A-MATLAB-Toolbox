function ITD = computeSphereITD(a,sPos,TYPE)
%COMPUTESPHEREITD Analytically computed ITD for a sphere.
%   B = COMPUTESPHEREITD(A,S) analytically computes ITD for a rigid sphere
%   of radius A (specified in meters), and for sound source positions, S, 
%   specified in SOFA cartesian coordinates, using the Woodworth formula 
%   (see, for example, Kuhn [1]). S must be specified as an N-by-3 matrix. 
%   B will contain ITD values in seconds.
%       If A is a scalar, B will be a row vector of length N.
%       If A is a vector of length M, B will be an M-by-N matrix.
%
%   B = COMPUTESPHEREITD(...,TYPE) optionally specifies the type of formula
%   to use to compute ITD. TYPE must be specified as a cell array. 
%   The options for TYPE are:
%       1. {'W'} - Woodworth formula (default). Antipodal ears are assumed.
%       2. {'LF'} - Low-frequency limit formula (see, for example, 
%       Kuhn [1]). Antipodal ears are assumed.
%       3. {'EW',[EL;ER],R} - Extended Woodworth formula. EL and ER
%       correspond to the left and right ear positions, respectively,
%       each specified as the coordinate pair (az,el) where az and el 
%       correspond to the azimuth and elevation, respectively, of the ear 
%       specified in degrees in SOFA spherical coordinates. R corresponds 
%       to source distance specified in meters. If R is not specified, R = 
%       inf is assumed. If EL and ER are not specified, EL = (100,-10) and 
%       ER = (260,-10) are assumed. The value of R specified here 
%       supersedes the source distance specified in S.

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

% Refs:
%   [1]. Kuhn (1977) - Model for the interaural time differences in the 
%   azimuthal plane.

narginchk(2,3);

if nargin < 3
    TYPE = {'W'};
end

% Check inputs
validateattributes(a,{'double'},{'vector','nonempty','nonnan','finite',...
    'positive','real'},'computeSphereITD','A',1);
validateattributes(sPos,{'double'},{'2d','nonempty','nonnan','finite',...
    'size',[NaN,3],'real'},'computeSphereITD','S',2);
validateattributes(TYPE,{'cell'},{'2d'},'computeSphereITD','TYPE',3);

% Get required variables and pre-process input data
a = shiftdim(a); % If a is a row vector, force it to be a column.
c = getSoundSpeed();

% Pre-process position data
switch lower(TYPE{1})
    case {'w','lf'}
        sDirs = sofaC2cipicI(sPos); % Convert to interaural coordinates
        theta = deg2rad(sDirs(:,1)).'; % Row vector of azimuths in radians
    case 'ew'
        if length(TYPE) < 3
            R = inf;
        else
            R = TYPE{3};
        end
        
        if length(TYPE) < 2
            EL = [100,-10];
            ER = [260,-10];
        else
            EL = TYPE{2}(1,:);
            ER = TYPE{2}(2,:);
        end
        
        [xEL,yEL,zEL] = sofaS2sofaC(EL(1),EL(2),1);
        [xER,yER,zER] = sofaS2sofaC(ER(1),ER(2),1);
        sDirs = sofaC2sofaS(sPos);
        S = sDirs(:,1:2);
        numPos = size(S,1);
        [xS,yS,zS] = sofaS2sofaC(S(:,1),S(:,2),ones(numPos,1));
        thetaL = getCentralAngle([xS,yS,zS],repmat([xEL,yEL,zEL],...
            numPos,1));
        thetaR = getCentralAngle([xS,yS,zS],repmat([xER,yER,zER],...
            numPos,1));
    otherwise
        error('Invalid TYPE{1} specification.')
end

% Main calculation
switch lower(TYPE{1})
    case 'w'
        ITD = a*(sin(theta)+theta)/c;
    case 'lf'
        ITD = 3*a*sin(theta)/c;
    case 'ew'
        rho = R/a;
        delL = zeros(numPos,1);
        delR = zeros(numPos,1);
        for ii = 1:numPos
            delL(ii) = getDel(thetaL(ii),rho);
            delR(ii) = getDel(thetaR(ii),rho);
        end
        ITD = (a/c)*(delL-delR);
    otherwise
        error('Invalid TYPE{1} specification.')
end

end

function D = getDel(T,R)
%GETDEL Compute delay value
%   D = GETDEL(T) computes the phase delay, D, for a given input incidence
%   angle, T, and non-dimensional source distance, R.

if R == inf
    if T <= 90
        D = -cosd(T);
    else
        D = deg2rad(T-90);
    end
else
    T0 = acosd(1/R);
    hatG = sqrt(R^2-1);
    G = sqrt(R^2-(2*R*cosd(T))+1);
    if T < T0
        D = G;
    else
        D = hatG + deg2rad(T-T0);
    end
end

end
