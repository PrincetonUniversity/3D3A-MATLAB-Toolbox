function varargout = sofaC2cipicI(varargin)
%SOFAC2CIPICI Convert SOFA cartesian to CIPIC interaural coordinates.
%   [A,E,R] = SOFAC2CIPICI(X,Y,Z) converts from SOFA cartesian coordinates, 
%   X, Y, and Z, to CIPIC interaural coordinates, A, E, and R. Azimuth, A, 
%   and elevation, E, are specified in degrees, while radius, R, is
%   specified in the same units as X, Y, and Z, which may each be specified
%   as scalars or vectors of the same length. A, E, and R will have the
%   same dimensions as X, Y, and Z.
%
%   CI = SOFAC2CIPICI(SC) allows X, Y, and Z to be specified as a 3-column 
%   matrix SC = [X,Y,Z].
%
%   ___ = SOFAC2CIPICI(___,'flipAz') optionally flips the sign convention
%   of A such that positions corresponding to positive values of Y in SOFA 
%   cartesian coordinates correspond to positive values of A (unlike in the 
%   CIPIC coordinate system, where these positions correspond to negative 
%   values).
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
            x = varargin{1};
            validateattributes(x,{'numeric'},{'2d','nonempty','finite',...
                'real','size',[NaN,3]},'sofaC2cipicI','SC',1)
            z = x(:,3);
            y = x(:,2);
            x = x(:,1);
        case 3
            singleArgFlag = false;
            x = varargin{1};
            y = varargin{2};
            z = varargin{3};
            
            validateattributes(x,{'numeric'},{'vector','finite','real'},...
                'sofaC2cipicI','X',1)
            validateattributes(y,{'numeric'},{'vector','finite','real',...
                'size',size(x)},'sofaC2cipicI','Y',2)
            validateattributes(z,{'numeric'},{'vector','finite','real',...
                'size',size(x)},'sofaC2cipicI','Z',3)
        otherwise
            error(['Invalid number/format of inputs given ''flipAz''',...
                ' is not specified.']) 
    end
    flipAzFlag = false;
else 
    switch indx
        case 2
            singleArgFlag = true;
            x = varargin{1};
            validateattributes(x,{'numeric'},{'2d','nonempty','finite',...
                'real','size',[NaN,3]},'sofaC2cipicI','SC',1)
            z = x(:,3);
            y = x(:,2);
            x = x(:,1);
        case 4
            singleArgFlag = false;
            x = varargin{1};
            y = varargin{2};
            z = varargin{3};
            
            validateattributes(x,{'numeric'},{'vector','finite','real'},...
                'sofaC2cipicI','X',1)
            validateattributes(y,{'numeric'},{'vector','finite','real',...
                'size',size(x)},'sofaC2cipicI','Y',2)
            validateattributes(z,{'numeric'},{'vector','finite','real',...
                'size',size(x)},'sofaC2cipicI','Z',3)
        otherwise
            error(['Invalid number/format of inputs given ''flipAz''',...
                ' is specified.'])
    end
    flipAzFlag = true;
end

rad = sqrt(x.^2 + y.^2 + z.^2);

if rad == 0
    error('X, Y, and Z cannot all be zero.')
end

el = mod(atan2d(z,x),360);
el(el >= 270) = el(el >= 270)-360;

if flipAzFlag
    az = -asind(-y./rad);
else
    az = asind(-y./rad);
end

if singleArgFlag
    if nargout > 1
        error('Only 1 output may be requested when input is SC.')
    else
        varargout{1} = [az,el,rad];
    end
else
    if nargout < 4
        varargout{1} = az;
        varargout{2} = el;
        varargout{3} = rad;
    else
        error('Too many output arguments requested.')
    end
end

end
