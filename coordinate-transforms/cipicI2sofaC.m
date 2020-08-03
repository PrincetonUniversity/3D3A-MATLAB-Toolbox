function varargout = cipicI2sofaC(varargin)
%CIPICI2SOFAC Convert CIPIC interaural to SOFA cartesian coordinates.
%   [X,Y,Z] = CIPICI2SOFAC(A,E,R) converts from CIPIC interaural 
%   coordinates, A (azimuth), E (elevation), and R (radius), to SOFA 
%   cartesian coordinates, X, Y, and Z. A, E, and R must be specified in 
%   degrees and may each be scalars or vectors of the same length. X, Y, Z
%   will have the same dimensions as A, E, and R, and will be specified in
%   the same units as R.
%
%   SC = CIPICI2SOFAC(CI) allows A, E, and R to be specified as a 3-column
%   matrix CI = [A,E,R]. X, Y, and Z are then returned as a 3-column
%   matrix, SC.
%
%   ___ = CIPICI2SOFAC(___,'flipAz') optionally indicates that the sign 
%   convention of A is flipped such that positions corresponding to 
%   positive values of Y in SOFA cartesian coordinates correspond to 
%   positive A (unlike in the CIPIC coordinate system, where these 
%   positions correspond to negative A).
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

narginchk(1,4);

indx = find(strcmpi(varargin,'flipAz'),1);
if isempty(indx)
    switch nargin
        case 1
            singleArgFlag = true;
            az = varargin{1};
            validateattributes(az,{'numeric'},{'2d','nonempty','finite',...
                'real','size',[NaN,3]},'cipicI2sofaC','CI',1)
            rad = az(:,3);
            el = az(:,2);
            az = az(:,1);
        case 3
            singleArgFlag = false;
            az = varargin{1};
            el = varargin{2};
            rad = varargin{3};
            
            validateattributes(az,{'numeric'},{'vector','finite',...
                'real'},'cipicI2sofaC','A',1)
            validateattributes(el,{'numeric'},{'vector','finite',...
                'real','size',size(az)},'cipicI2sofaC','E',2)
            validateattributes(rad,{'numeric'},{'vector','finite',...
                'real','size',size(az)},'cipicI2sofaC','R',3)
        otherwise
            error(['Invalid number/format of inputs given ''flipAz''',...
                ' is not specified.']) 
    end
    flipAzFlag = false;
else 
    switch indx
        case 2
            singleArgFlag = true;
            az = varargin{1};
            validateattributes(az,{'numeric'},{'2d','nonempty','finite',...
                'real','size',[NaN,3]},'cipicI2sofaC','CI',1)
            rad = az(:,3);
            el = az(:,2);
            az = az(:,1);
        case 4
            singleArgFlag = false;
            az = varargin{1};
            el = varargin{2};
            rad = varargin{3};
            
            validateattributes(az,{'numeric'},{'vector','finite',...
                'real'},'cipicI2sofaC','A',1)
            validateattributes(el,{'numeric'},{'vector','finite','real',...
                'size',size(az)},'cipicI2sofaC','E',2)
            validateattributes(rad,{'numeric'},{'vector','finite',...
                'real','size',size(az)},'cipicI2sofaC','R',3)
        otherwise
            error(['Invalid number/format of inputs given ''flipAz''',...
                ' is specified.'])
    end
    flipAzFlag = true;
end

validateattributes(az,{'numeric'},{'vector','>=',-90,'<=',90},...
    'cipicI2sofaC','A')
validateattributes(el,{'numeric'},{'vector','>=',-90,'<',270},...
    'cipicI2sofaC','E')
validateattributes(rad,{'numeric'},{'vector','positive'},'cipicI2sofaC',...
    'R')

if flipAzFlag
    az = -az;
end

x = rad.*cosd(-az).*cosd(el);
y = rad.*sind(-az);
z = rad.*cosd(-az).*sind(el);

if singleArgFlag
    if nargout > 1
        error('Only 1 output may be requested when input is CI.')
    else
        varargout{1} = [x,y,z];
    end
else
    if nargout < 4
        varargout{1} = x;
        varargout{2} = y;
        varargout{3} = z;
    else
        error('Too many output arguments requested.')
    end
end

end
