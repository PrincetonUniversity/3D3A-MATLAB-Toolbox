function varargout = auditoryParallax(pIn,ePos,varargin)
%AUDITORYPARALLAX Auditory parallax estimation.
%   PO = AUDITORYPARALLAX(PI,E) returns output positions, PO, with
%   reference to a coordinate system centered around E given input
%   positions positions, PI, in a coordinate centered around the origin.
%   All positions must be specified in SOFA cartesian coordinates. E must
%   be specified with respect to the same origin as PI. PI must be
%   specified as an N-by-3 matrix, where N is the number of positions, 
%   while E must be specified as a 1-by-3 vector. PO will have the same 
%   dimensions as PI.
%
%   [PO,PIH] = AUDITORYPARALLAX(PI,E,R) optionally specifies the radius, R,
%   at which the parallax should be computed and PIH contains the
%   positions, PI, when projected onto a circle with radius R centered at
%   the same origin as PI by extending a straight line from E through PI to
%   R. For more, see the "cross-ear" method described by Romblom and Cook
%   [1]. The default value for R is inf, for which PIH is equal to PO.
%
%   [PO,PIH] = AUDITORYPARALLAX(PI,E,R,METHOD) optionally specifies the
%   method to use. METHOD can take the following options:
%       1. 'ear' (default) - standard approach
%       2. 'tangent' - trial approach

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

% References:
%   [1]. Romblom and Cook (2008) - Near-Field Compensation for HRTF 
%   Processing.

narginchk(2,4);

% Check inputs
validateattributes(pIn,{'numeric'},{'2d','nonempty','nonnan','finite',...
    'real','size',[NaN,3]},'auditoryParallax','PI',1);
validateattributes(ePos,{'numeric'},{'2d','nonempty','nonnan','finite',...
    'size',[1,3],'real'},'auditoryParallax','E',2);

if nargin < 4
    method = 'ear';
else
    method = varargin{2};
    validateattributes(method,{'char'},{'scalartext'},...
        'auditoryParallax','method',4);
end

if nargin < 3
    pihFlag = false;
else
    R = varargin{1};
    if R == inf || isempty(R)
        pihFlag = false;
    else
        validateattributes(R,{'numeric'},{'scalar','nonnan','finite',...
            'positive'},'auditoryParallax','R',3);
        pihFlag = true;
    end
end

numPos = size(pIn,1);
switch lower(method)
    case 'ear'
        pO = pIn-repmat(ePos,numPos,1);
    case 'tangent'
        eDirs = sofaC2sofaS(ePos);
        pIn_dirs = sofaC2sofaS(pIn);
        rho = pIn_dirs(1,3)/eDirs(1,3);
        ePosMat = repmat(ePos,numPos,1);
        theta = getCentralAngle(pIn,ePosMat);
        theta_sub1 = theta <= acosd(1/rho);
        theta_sub2 = theta > acosd(1/rho);
        sub1_count = length(find(theta_sub1));
        sub2_count = length(find(theta_sub2));
        
        pO = pIn; % Initialization
        pO(theta_sub1,:) = pIn(theta_sub1,:)-repmat(ePos,sub1_count,1);
        
        vert_angle = acosd(eDirs(1,3)/pIn_dirs(1,3));
        pIn_int_dirs = sofaC2cipicI(pIn,'sofaAz'); % interaural coordinates
        angle_opt_1 = mod(pIn_int_dirs(theta_sub2,1)-vert_angle,360);
        angle_opt_2 = mod(pIn_int_dirs(theta_sub2,1)+vert_angle,360);
        rad_int_dirs_1 = [angle_opt_1,pIn_int_dirs(theta_sub2,2),...
            repmat(eDirs(1,3),sub2_count,1)];
        rad_int_dirs_2 = [angle_opt_2,pIn_int_dirs(theta_sub2,2),...
            repmat(eDirs(1,3),sub2_count,1)];
        rad_pos_1 = cipicI2sofaC(rad_int_dirs_1,'sofaAz');
        rad_pos_2 = cipicI2sofaC(rad_int_dirs_2,'sofaAz');
        
        rad_pos = rad_pos_1; % Initialize
        check_angle_1 = getCentralAngle(rad_pos_1,ePosMat(theta_sub2,:));
        check_angle_2 = getCentralAngle(rad_pos_2,ePosMat(theta_sub2,:));
        indxList = check_angle_1 <= check_angle_2;
        rad_pos(indxList,:) = rad_pos_1(indxList,:);
        indxList = check_angle_1 > check_angle_2;
        rad_pos(indxList,:) = rad_pos_2(indxList,:);
        pO(theta_sub2,:) = pIn(theta_sub2,:)-rad_pos;
    otherwise
        error('Invalid method specification.')
end

if pihFlag
    pO_dirs = sofaC2sofaS(pO);
    ePosMat = repmat(ePos,numPos,1);
    eDirMat = sofaC2sofaS(ePosMat);
    theta = getCentralAngle(pIn,ePosMat);
    rE = sqrt(R.^2-(2.*eDirMat(:,3).*R.*cosd(theta))+(eDirMat(:,3).^2));
    pO_dirs(:,3) = rE;
    pO_s = sofaS2sofaC(pO_dirs);
    pIH = ePosMat+pO_s;
else
    pIH = pO;
end

switch nargout
    case 1
        varargout{1} = pO;
    case 2
        varargout{1} = pO;
        varargout{2} = pIH;
    otherwise
        error('Too many outputs requested.')
end

end
