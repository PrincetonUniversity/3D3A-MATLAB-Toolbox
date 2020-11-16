function varargout = sofaS2sofaC(varargin)
%SOFAS2SOFAC Convert SOFA spherical to SOFA cartesian coordinates.
%   [X,Y,Z] = SOFAS2SOFAC(A,E,R) converts from SOFA spherical to SOFA
%   cartesian coordinates. Input azimuth (A) and elevation (E) must be 
%   specified in degrees, where A can be in the range [0,360) or
%   (-180,180], and E in the range [-90,90]. A, E, and R may be specified 
%   as scalars or vectors. If vectors, they must have the same length. X, 
%   Y, and Z will have the same units as R and the same dimensions as A, E, 
%   and R.
%
%   SC = SOFAS2SOFAC(SS) allows A, E, and R to be specified as a 3-column 
%   matrix [A,E,R]. SC is then the 3-column matrix [X,Y,Z].
%
%   See also SOFAC2SOFAS, CIPICI2SOFAC, SOFAC2CIPICI.

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

narginchk(1,3)

switch nargin
    case 1
        singleArgFlag = true;
        az = varargin{1};
        validateattributes(az,{'numeric'},{'2d','nonempty','finite',...
            'real','size',[NaN,3]},'sofaS2sofaC','SS',1)
        rad = az(:,3);
        el = az(:,2);
        az = az(:,1);
    case 3
        singleArgFlag = false;
        az = varargin{1};
        el = varargin{2};
        rad = varargin{3};
        
        validateattributes(az,{'numeric'},{'vector','finite','real'},...
            'sofaS2sofaC','A',1)
    otherwise
        error('Invalid number/format of inputs.')
end

validateattributes(az,{'numeric'},{'vector','>=',-180,'<',360},...
    'sofaS2sofaC','A')
validateattributes(el,{'numeric'},{'vector','finite','real',...
    'size',size(az),'>=',-90,'<=',90},'sofaS2sofaC','E')
validateattributes(rad,{'numeric'},{'vector','finite',...
    'real','size',size(az),'positive'},'sofaS2sofaC','R')

x = rad.*cosd(el).*cosd(az);
y = rad.*cosd(el).*sind(az);
z = rad.*sind(el);

if singleArgFlag
    if nargout > 1
        error('Only 1 output may be requested when input is SS.')
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
