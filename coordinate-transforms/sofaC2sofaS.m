function varargout = sofaC2sofaS(varargin)
%SOFAC2SOFAS Convert SOFA cartesian to SOFA spherical coordinates.
%   [A,E,R] = SOFAC2SOFAS(X,Y,Z) converts from SOFA cartesian to SOFA
%   spherical coordinates. X, Y, and Z may be scalars or vectors with the 
%   same dimensions. A (azimuth), E (elevation), and R (radius) will have 
%   the same dimensions as X, Y, and Z. A and E will be specified in 
%   degrees where A will be in the range [0,360) and E in the range 
%   [-90,90]. R will have the same units as X, Y, and Z.
%
%   SS = SOFAC2SOFAS(SC) allows X, Y, and Z to be specified as a 3-column
%   matrix [X,Y,Z]. SS is then the 3-column matrix [A,E,R].
%
%   ___ = SOFAC2SOFAS(___,'azRange',RANGE) optionally specifies the range 
%   for A. The two possible options are:
%       1. 'pos' (default) - A ranges from 0 to 360 (not inclusive).
%
%       2. 'pm' - A ranges from -180 (not inclusive) to 180.
%
%   See also SOFAS2SOFAC, CIPICI2SOFAC, SOFAC2CIPICI.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   Joseph G. Tylka <josephgt@princeton.edu>
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

narginchk(1,5)

indx = find(strcmpi(varargin,'azRange'),1);
if isempty(indx)
    switch nargin
        case 1
            singleArgFlag = true;
            x = varargin{1};
            validateattributes(x,{'numeric'},{'2d','nonempty','finite',...
                'real','size',[NaN,3]},'sofaC2sofaS','SC',1)
            z = x(:,3);
            y = x(:,2);
            x = x(:,1);
        case 3
            singleArgFlag = false;
            x = varargin{1};
            y = varargin{2};
            z = varargin{3};
            
            validateattributes(x,{'numeric'},{'vector','finite','real'},...
                'sofaC2sofaS','X',1)
            validateattributes(y,{'numeric'},{'vector','finite','real',...
                'size',size(x)},'sofaC2sofaS','Y',2)
            validateattributes(z,{'numeric'},{'vector','finite','real',...
                'size',size(x)},'sofaC2sofaS','Z',3)
        otherwise
            error(['Invalid number/format of inputs given ''azRange''',...
                ' is not specified.']) 
    end
    RANGE = 'pos'; % Default value for R when 'azRange' is not specified.
else 
    if length(varargin) ~= (indx+1)
        error('''azRange'' option must be followed by a value for RANGE.');
    else
        RANGE = varargin{indx+1};
        validateattributes(RANGE,{'char'},{'scalartext','nonempty'},...
            'sofaC2sofaS','RANGE')
    end
    
    switch indx
        case 2
            singleArgFlag = true;
            x = varargin{1};
            validateattributes(x,{'numeric'},{'2d','nonempty','finite',...
                'real','size',[NaN,3]},'sofaC2sofaS','SC',1)
            z = x(:,3);
            y = x(:,2);
            x = x(:,1);
        case 4
            singleArgFlag = false;
            x = varargin{1};
            y = varargin{2};
            z = varargin{3};
            
            validateattributes(x,{'numeric'},{'vector','finite','real'},...
                'sofaC2sofaS','X',1)
            validateattributes(y,{'numeric'},{'vector','finite','real',...
                'size',size(x)},'sofaC2sofaS','Y',2)
            validateattributes(z,{'numeric'},{'vector','finite','real',...
                'size',size(x)},'sofaC2sofaS','Z',3)
        otherwise
            error(['Invalid number/format of inputs given ''azRange''',...
                ' is specified.'])
    end
end

rad = sqrt(x.^2 + y.^2 + z.^2);

if rad == 0
    error('X, Y, and Z cannot all be zero.')
end

switch lower(RANGE)
    case 'pos'
        % Two mods are required because mod(-1e-14,360) returns 360 instead
        % of 0.
        az = mod(mod(atan2d(y,x),360),360);
    case 'pm'
        az = atan2d(y,x);
    otherwise
        error('Invalid RANGE specification.')
end

el = asind(z./rad);

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
