function ITD = computeSphereITD(a,sPos,TYPE)
%COMPUTESPHEREITD Analytically computed ITD for a sphere.
%   B = COMPUTESPHEREITD(A,S) analytically computes ITD for a rigid sphere
%   of radius A (specified in meters), and for sound source positions, S, 
%   specified in SOFA cartesian coordinates, using the Woodworth and 
%   Schlosberg formula (see, for example, Kuhn [1]). S must be specified as 
%   an N-by-3 matrix. B will contain ITD values in seconds.
%       If A is a scalar, B will be a row vector of length N.
%       If A is a vector of length M, B will be an M-by-N matrix.
%
%   B = COMPUTESPHEREITD(...,TYPE) optionally specifies the type of formula
%   to use to compute ITD. The two options for TYPE are:
%       1. 'WS' - Woodworth and Schlosberg formula (default).
%       2. 'LF' - Low-frequency limit formula (see, for example, Kuhn [1]).

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

% Refs:
%   [1]. Kuhn (1977) - Model for the interaural time differences in the 
%   azimuthal plane.

narginchk(2,3);

if nargin < 3
    TYPE = 'WS';
end

% Check inputs
validateattributes(a,{'double'},{'vector','nonempty','nonnan','finite',...
    'positive','real'},'computeSphereITD','A',1);
validateattributes(sPos,{'double'},{'2d','nonempty','nonnan','finite',...
    'size',[NaN,3],'real'},'computeSphereITD','S',2);
validateattributes(TYPE,{'char'},{'scalartext','nonempty'},...
    'computeSphereITD','TYPE',3);

% Reformat inputs
a = shiftdim(a); % If a is a row vector, force it to be a column.
sDir = sofaC2cipicI(sPos); % Convert to interaural coordinates
theta = deg2rad(sDir(:,1)).'; % Row vector of azimuths in radians

% Main calculation
c = getSoundSpeed();
switch lower(TYPE)
    case 'ws'
        ITD = a*(sin(theta)+theta)/c;
    case 'lf'
        ITD = 3*a*sin(theta)/c;
    otherwise
        error('Invalid TYPE specification.')
end

end
