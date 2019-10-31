function varargout = sofaC2cipicI(varargin)
%SOFAC2CIPICI Convert SOFA cartesian coordinates to CIPIC interaural 
%coordinates.
%   [A,E,R] = SOFAC2CIPICI(X,Y,Z) converts from SOFA cartesian coordinates, 
%   X, Y, and Z, to CIPIC interaural coordinates, A, E, and R. Azimuth, A, 
%   and elevation, E, are specified in degrees, while radius, R, is
%   specified in the same units as X, Y, and Z.
%       If X, Y, and Z are specified as scalars, A, E, and R will also be
%       scalars.
%       If X, Y, and X are specified as vectors, they must all have the
%       same length, N. A, E, and R will then be column vectors of length
%       N.
%
%   CI = SOFAC2CIPICI(X,Y,Z) optionally returns an N-by-3 matrix, CI = 
%   [A,E,R].
%
%   __ = SOFAC2CIPICI(SC) allows X, Y, and Z to be specified as a 3-column 
%   matrix SC = [X,Y,Z].
%
%   __ = SOFAC2CIPICI(__,'flipAz') optionally flips the sign convention of
%   azimuth values such that positions corresponding to positive values of
%   Y in SOFA cartesian coordinates correspond to positive azimuths
%   (unlike in the CIPIC coordinate system, where these positions
%   correspond to negative azimuths).
%
%   __ = SOFAC2CIPICI(__,'sofaAz') optionally returns azimuths in SOFA
%   spherical coordinates ranging from 0 to 360 degrees. Note that the 
%   'flipAz' option has no effect if this option is also specified.
%
%   See also CIPICI2SOFAC, SOFAS2SOFAC, SOFAC2SOFAS.

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

narginchk(1,5);

validateattributes(varargin{1},{'double'},{'2d','nonempty','nonnan',...
    'finite','real'},'sofaC2cipicI','X or SC',1);
switch size(varargin{1},2)
    case 1
        x = shiftdim(varargin{1});
        validateattributes(x,{'double'},{'vector','nonempty','nonnan',...
            'finite','real'},'sofaC2cipicI','X',1);
        xLen = length(x);
        y = shiftdim(varargin{2});
        validateattributes(y,{'double'},{'vector','nonempty','nonnan',...
            'finite','real','numel',xLen},'sofaC2cipicI','Y',2);
        z = shiftdim(varargin{3});
        validateattributes(z,{'double'},{'vector','nonempty','nonnan',...
            'finite','real','numel',xLen},'sofaC2cipicI','Z',3);
    case 3
        SC = varargin{1};
        validateattributes(SC,{'double'},{'2d','nonempty','nonnan',...
            'finite','real','size',[NaN,3]},'sofaC2cipicI','SC',1);
        x = SC(:,1);
        y = SC(:,2);
        z = SC(:,3);
    otherwise
        error('First input has invalid dimensions.')
end

rad = sqrt(x.^2 + y.^2 + z.^2);

if rad == 0
    error('X, Y, and Z cannot all be zero.')
end

el = mod(atan2d(z,x),360);
el(el >= 270) = el(el >= 270)-360;

flag1 = find(strcmpi(varargin,'sofaAz'),1);
az = el; % Initialize
if flag1
    subIndxs = el > 90;
    az(subIndxs) = 180+asind(-y(subIndxs)./rad(subIndxs));
    subIndxs = el <= 90;
    az(subIndxs) = mod(-asind(-y(subIndxs)./rad(subIndxs)),360);
else
    flag2 = find(strcmpi(varargin,'flipAz'),1);
    if flag2
        az = -asind(-y./rad);
    else
        az = asind(-y./rad);
    end
end

switch nargout
    case {0,1}
        varargout{1} = [az,el,rad];
    case 2
        varargout{1} = az;
        varargout{2} = el;
    case 3
        varargout{1} = az;
        varargout{2} = el;
        varargout{3} = rad;
    otherwise
        error('Invalid number of requested outputs.')
end

end
