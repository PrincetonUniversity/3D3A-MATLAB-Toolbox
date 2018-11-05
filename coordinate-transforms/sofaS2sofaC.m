function [x,y,z] = sofaS2sofaC(az,el,rad)
%SOFAS2SOFAC Convert SOFA spherical coordinates to SOFA cartesian
%coordinates.
%   [x,y,z] = SOFAS2SOFAC(az,el,rad) converts from SOFA spherical 
%   coordinates to SOFA cartesian coordinates. Input azimuth (az) and 
%   elevation (el) must be specified in degrees. az, el, and rad may be 
%   specified as scalars or vectors. If vectors, they must have the same 
%   length. x, y, and z are either scalars, or column vectors with the same
%   length as az, el, and rad.
%
%   SC = SOFAS2SOFAC(SS) allows az, el, and rad to be specified as a 
%   3-column matrix [az,el,rad]. SC is then the 3-column matrix [x,y,z].
%
%   See also SOFAC2SOFAS.

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

narginchk(1,3);

singleArgFlag = false;
if nargin == 1 && nargout <= 1
    singleArgFlag = true;
    
    validateattributes(az,{'double'},{'2d','nonempty','nonnan','finite',...
        'real','size',[NaN,3]},'sofaS2sofaC','SS',1)
    
    rad = az(:,3);
    el = az(:,2);
    az = az(:,1);
elseif nargin == 3
    
    validateattributes(az,{'double'},{'vector','nonempty','nonnan',...
        'finite','real'},'sofaS2sofaC','az',1)
    validateattributes(el,{'double'},{'vector','nonempty','nonnan',...
        'finite','real','numel',numel(az)},'sofaS2sofaC','el',2)
    validateattributes(rad,{'double'},{'vector','nonempty','nonnan',...
        'finite','real','numel',numel(az)},'sofaS2sofaC','rad',3)
    
    az = shiftdim(az);
    el = shiftdim(el);
    rad = shiftdim(rad);
else
    error('2 inputs provided when either 1 or 3 inputs are required.')
end

if (max(az) >= 360) || (min(az) < 0)
    error('One or more az values are outside the range [0,360).')
end

if (max(el) > 90) || (min(el) < -90)
    error('One or more el values are outside the range [-90,90].')
end

x = rad.*cosd(el).*cosd(az);
y = rad.*cosd(el).*sind(az);
z = rad.*sind(el);

if singleArgFlag
    x = [x y z];
end

end
