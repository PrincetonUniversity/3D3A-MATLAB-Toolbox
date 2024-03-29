function psi = getArrayPotential(a,k,r)
%GETARRAYPOTENTIAL Potential due to an array of sources.
%   PSI = GETARRAYPOTENTIAL(A,K,R) computes the potential field at points,
%   R, in space and for a given angular wavenumber, K, generated by an
%   array of sound sources whose properties are specified in the structure,
%   A. K must be specified in rad/m and must be a vector of length N, where 
%   N >= 1. R must be specified in meters using Cartesian coordinates, and 
%   must be an M-by-3 matrix where M >= 1.
%
%   The following fields need to be specified in the structure, A:
%       1. 'sourceType' - a character array that can have the following
%       values:
%           (i) 'point source' - each source on the array is a point
%           source.
%           (ii) 'baffled circular piston' (default) - each source on the 
%           array is a baffled circular piston which is simulated by 
%           assuming that the potential field of interest is in the far-
%           field of each source (i.e., the potential field is computed at 
%           a distance that is at least an order of magnitude larger than 
%           the piston radius).
%
%       2. 'numSources' - the number of sources on the array. This must be
%       a positive integer.
%
%       3. 'sourceRadii' - the radius of each source if sourceType is
%       'baffled circular piston'. This field is ignored if sourceType is
%       'point source'. This must be a numSources-length vector.
%
%       4. 'sourcePos' - the positions of the sources on the array in
%       meters in global Cartesian coordinates. This must be a numSources-
%       by-3 matrix.
%
%       5. 'sourceAxes' - the axes of each source on the array in global 
%       Cartesian coordinates.
%
%       6. 'pistonAccel' - scalar piston acceleration in m/s^2, required if
%       sourceType is 'baffled circular piston'. The default value is 
%       1/(1.21*(sourceRadii)^2). This must be a numSources-length vector.
%
%       7. 'addDelay' - (optional) additional delay in seconds. The default
%       value is 0. This must be a numSources-length vector.
%
%   The output potential field, PSI, will have dimensions of N-by-M. The 
%   GETPRESSURE function may be used to compute the pressure from the 
%   returned PSI value(s).
%
%   See also GETPOINTSOURCEPOTENTIAL, GETPISTONPOTENTIAL, GETPRESSURE.

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

% Check number of inputs
narginchk(3,3);

validateattributes(a,{'struct'},{},'getArrayPotential','P',1);
validateattributes(k,{'numeric'},{'vector','real','finite',...
    'nonnegative'},'getArrayPotential','K',2);
validateattributes(r,{'numeric'},{'2d','real','finite','size',[NaN,3]},...
    'getArrayPotential','R',3);

k = shiftdim(k).'; % Force k to be a row vector
kLen = length(k);
numPos = size(r,1);

if ~isfield(a,'sourceType')
    a.sourceType = 'baffled circular piston';
end
validateattributes(a.sourceType,{'char'},{'scalartext'},...
    'getArrayPotential','A.SOURCETYPE',1);

if ~isfield(a,'numSources')
    error('Required field ''numSources'' not present in A.')
end
validateattributes(a.numSources,{'numeric'},{'scalar','real','finite',...
    'positive'},'getArrayPotential','A.NUMSOURCES',1);

if ~isfield(a,'sourceRadii')
    error('Required field ''sourceRadii'' not present in A.')
end
validateattributes(a.sourceRadii,{'numeric'},{'vector','real','finite',...
    'numel',a.numSources},'getArrayPotential','A.SOURCERADII',1);
a.sourceRadii = shiftdim(a.sourceRadii);

if ~isfield(a,'sourcePos')
    error('Required field ''sourcePos'' not present in A.')
end
validateattributes(a.sourcePos,{'numeric'},{'2d','real','finite',...
    'size',[a.numSources,3]},'getArrayPotential','A.SOURCEPOS',1);

if ~isfield(a,'sourceAxes')
    error('Required field ''sourceAxes'' not present in A.')
end
validateattributes(a.sourceAxes,{'numeric'},{'2d','real','finite',...
    'size',[a.numSources,3]},'getArrayPotential','A.SOURCEAXES',1);

if ~isfield(a,'pistonAccel')
    a.pistonAccel = 1./(1.21*(a.sourceRadii).^2);
end
validateattributes(a.pistonAccel,{'numeric'},{'vector','real','finite',...
    'numel',a.numSources},'getArrayPotential','A.PISTONACCEL',1);
a.pistonAccel = shiftdim(a.pistonAccel);

if ~isfield(a,'addDelay')
    a.addDelay = zeros(a.numSources,1);
end
validateattributes(a.addDelay,{'numeric'},{'vector','real','finite',...
    'numel',a.numSources},'getArrayPotential','A.ADDDELAY',1);
a.addDelay = shiftdim(a.addDelay);

psi = zeros(kLen,numPos);
switch lower(a.sourceType)
    case 'point source'
        for ii = 1:(a.numSources)
            piston.position = a.sourcePos(ii,:);
            piston.radius = a.sourceRadii(ii);
            piston.axis = a.sourceAxes(ii,:);
            currPsi = getPistonPotential(a.sourcePos(ii,:),k,r,...
                a.addDelay(ii));
            psi = psi + currPsi;
        end
    case 'baffled circular piston'
        for ii = 1:(a.numSources)
            piston.position = a.sourcePos(ii,:);
            piston.radius = a.sourceRadii(ii);
            piston.axis = a.sourceAxes(ii,:);
            currPsi = getPistonPotential(piston,k,r,a.pistonAccel(ii),...
                a.addDelay(ii));
            psi = psi + currPsi;
        end
    otherwise
        error('Invalid sourceType specification.')
end

end
