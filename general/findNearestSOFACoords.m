function [indx,err] = findNearestSOFACoords(posMat,target,COORDSYS,DIST)
%FINDNEARESTSOFACOORDS Find the nearest specified coordinates from given 
%set of SOFA coordinates and return the corresponding index.
%   [indx,err] = FINDNEARESTSOFACOORDS(posMat,target) finds the coordinates
%   in posMat to which target is nearest (in an l2 sense) and returns the
%   corresponding index along with the error. posMat must be an N-by-3 
%   matrix specifying positions in SOFA cartesian coordinates. target must 
%   be a vector of length 3 specifying the target position in SOFA 
%   cartesian coordinates.
%
%   [indx,err] = FINDNEARESTSOFACOORDS(...,COORDSYS) optionally specifies 
%   the coordinate system used to specify target. posMat must always be in 
%   SOFA cartesian coordinates. Options for COORDSYS are:
%       1. 'SOFAc' (default) - SOFA cartesian coordinates
%       2. 'SOFAs' - SOFA spherical coordinates
%       3. 'CIPICi' - CIPIC interaural coordinates
%
%   [indx,err] = FINDNEARESTSOFACOORDS(...,DISTSPEC) optionally specifies 
%   whether or not the distance from the origin should be considered when
%   determining the nearest point. The two options are: 'on' and 'off'
%   (default).

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
%   Copyright (c) 2018 Princeton University
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

if nargin < 4
    DIST = 'off';
end

if nargin < 3
    COORDSYS = 'SOFAc';
end

[numRows,numCols] = size(posMat);
if numCols ~= length(target) || numCols ~= 3
    error('posMat must be N-by-3 and target 1-by-3')
else
    target = reshape(target,[1,3]);
end

switch lower(COORDSYS)
    case 'sofac'
        targetC = target;
    case 'sofas'
        targetC = sofaS2sofaC(target);
    case 'cipici'
        [x,y,z] = cipic2sofaC(target(1),target(2),target(3));
        targetC = [x,y,z];
    otherwise
        error('Invalid specification for COORDSYS.')
end

if strcmpi(DIST,'off')
    normMat = sqrt(sum(posMat.^2,2))*ones(1,numCols);
    posMat = posMat./normMat;
    targetC = targetC/sqrt(sum(targetC.^2));
end

distVec = sqrt(sum((posMat-ones(numRows,1)*targetC).^2,2));
[err,indx] = min(distVec);

end

