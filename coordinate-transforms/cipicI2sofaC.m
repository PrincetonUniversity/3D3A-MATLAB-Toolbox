function varargout = cipicI2sofaC(varargin)
%CIPICI2SOFAC Convert CIPIC interaural coordinates to SOFA cartesian
%coordinates.
%   [X,Y,Z] = CIPICI2SOFAC(A,E,R) converts from CIPIC interaural 
%   coordinates, A, E, and R, to SOFA cartesian coordinates, X, Y, and Z. 
%   Input azimuth, A, and elevation, E, must be specified in degrees. 
%       If A, E, and R are scalars, then X, Y, and Z are also scalars
%       specified in the same units as R.
%       If A, E, and R, are vectors, they must all have the same length, N.
%       Then, X, Y, and Z will be column vectors of length, N.
%
%   SC = CIPICI2SOFAC(A,E,R) optionally returns an N-by-3 matrix, SC = 
%   [X,Y,Z].
%
%   __ = CIPICI2SOFAC(CI) allows A, E, and R to be specified as a 3-column
%   matrix [A,E,R].
%
%   __ = CIPICI2SOFAC(__,'flipAz') optionally indicates that the sign 
%   convention of azimuth is flipped such that positions corresponding to 
%   positive values of Y in SOFA cartesian coordinates correspond to 
%   positive azimuths (unlike in the CIPIC coordinate system, where these 
%   positions correspond to negative azimuths).
%
%   __ = CIPICI2SOFAC(__,'sofaAz') optionally allows azimuths to be
%   specified in SOFA spherical coordinates ranging from 0 to 360 degrees. 
%   Note that the 'flipAz' option has no effect if this option is also 
%   specified.
%
%   See also SOFAC2CIPICI, SOFAC2SOFAS, SOFAS2SOFAC.

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
    'finite','real'},'cipicI2sofaC','A or CI',1);
switch size(varargin{1},2)
    case 1
        az = shiftdim(varargin{1});
        validateattributes(az,{'double'},{'vector','nonempty','nonnan',...
            'finite','real'},'cipicI2sofaC','A',1);
        azLen = length(az);
        el = shiftdim(varargin{2});
        validateattributes(el,{'double'},{'vector','nonempty','nonnan',...
            'finite','real','numel',azLen},'cipicI2sofaC','E',2);
        rad = shiftdim(varargin{3});
        validateattributes(rad,{'double'},{'vector','nonempty','nonnan',...
            'finite','real','numel',azLen},'cipicI2sofaC','R',3);
    case 3
        CI = varargin{1};
        validateattributes(CI,{'double'},{'2d','nonempty','nonnan',...
            'finite','real','size',[NaN,3]},'cipicI2sofaC','CI',1);
        az = CI(:,1);
        el = CI(:,2);
        rad = CI(:,3);
    otherwise
        error('First input has invalid dimensions.')
end

if (max(el) >= 270) || (min(el) < -90)
    error('One or more E values are outside the range [-90,270).')
end

flag1 = find(strcmpi(varargin,'sofaAz'),1);
if flag1
    if (min(az) < 0) || (max(az) >= 360)
        error('One or more az values are outside the range [0,360).')
    end
    
    subIndxs = (el >= -90) & (el <= 90) & (az >= 270);
    az(subIndxs) = az(subIndxs)-360;
    subIndxs = (el > 90) & (el < 270);
    az(subIndxs) = 180-az(subIndxs);
else
    if (max(az) > 90) || (min(az) < -90)
        error('One or more az values are outside the range [-90,90].')
    end
    
    flag2 = find(strcmpi(varargin,'flipAz'),1);
    if flag2
        az = -az;
    end
end
x = rad.*cosd(-az).*cosd(el);
y = rad.*sind(-az);
z = rad.*cosd(-az).*sind(el);

switch nargout
    case {0,1}
        varargout{1} = [x,y,z];
    case 2
        varargout{1} = x;
        varargout{2} = y;
    case 3
        varargout{1} = x;
        varargout{2} = y;
        varargout{3} = z;
    otherwise
        error('Invalid number of requested outputs.')
end

end
