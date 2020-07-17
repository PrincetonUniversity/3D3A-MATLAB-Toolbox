function ILD = computeSphereILD(a,sPos,varargin)
%COMPUTESPHEREILD Exact low-frequency ILD using a rigid-sphere model.
%   B = COMPUTESPHEREILD(A,S) analytically computes low-frequency ILD for a 
%   rigid sphere of radius A (specified in meters), and for sound source 
%   positions, S, specified in SOFA cartesian coordinates, using the 
%   exact formula derived by Sridhar and Choueiri. S must be specified 
%   as an N-by-3 matrix. B will contain ILD values in dB. It is assumed
%   that the ears are antipodal.
%       If A is a scalar, B will be a row vector of length N.
%       If A is a vector of length M, B will be an M-by-N matrix.
%
%   B = COMPUTESPHEREILD(...,[EL;ER]) optionally specifies the positions of
%   the left and right "ears" on the surface of the sphere. EL and ER
%   correspond to the left and right ear positions, respectively, each
%   specified as the coordinate pair (az,el) where az and el correspond to
%   the azimuth and elevation, respectively, of the ear specified in 
%   degrees in SOFA spherical coordinates. 
%
%   B = COMPUTESPHEREILD(...,[EL;ER],R) optionally specifies R as the 
%   source distance in meters. If R is not specified, R is computed 
%   directly from S. If R is specified, its value supersedes the source
%   distance corresponding to the source positions given in S. R must be
%   specified in meters. For a source infinitely far away, specify R as
%   inf.

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
    'positive','real'},'computeSphereILD','A',1);
validateattributes(sPos,{'numeric'},{'2d','nonempty','nonnan','finite',...
    'size',[NaN,3],'real'},'computeSphereILD','S',2);

% Get required variables and pre-process input data
a = shiftdim(a); % If a is a row vector, force it to be a column.
sDirs = sofaC2sofaS(sPos);
numDirs = size(sDirs,1); % numDirs = N in documentation.

if nargin < 4
    R = sDirs(:,3);
else
    R = varargin{2};
    validateattributes(R,{'numeric'},{'vector','nonempty','nonnan'},...
        'computeSphereILD','R',4);
    
    if isscalar(R)
        R = R*ones(1,numDirs);
    else
        if length(R) ~= numDirs
            error('Invalid R specification.')
        end
    end
end

if nargin < 3
    EL = [90,0];
    ER = [270,0];
else
    ePos = varargin{1};
    validateattributes(ePos,{'numeric'},{'2d','size',[2,2],'real'},...
        'computeSphereILD','[EL;ER]',3);
    
    EL = ePos(1,:);
    ER = ePos(2,:);
end

% DC gain calculation
HL = computeSphereLFGain(a,sPos,EL,R);
HR = computeSphereLFGain(a,sPos,ER,R);

% Ratio of DC gains, in dB
ILD = mag2db(HL./HR);

end
