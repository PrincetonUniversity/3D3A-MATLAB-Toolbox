function [x,y,z] = cipic2sofaC(az,el,rad)
%CIPIC2SOFAC Convert CIPIC interaural polar coordinates to SOFA cartesian
%coordinates.
%   [x,y,z] = CIPIC2SOFAC(az,el,rad) converts from CIPIC interaural polar 
%   coordinates to SOFA cartesian coordinates. Input azimuth (az) and 
%   elevation (el) must be specified in degrees. az, el, and rad may be 
%   specified as scalars or vectors. If vectors, they must have the same 
%   length. x, y, and z are either scalars or column vectors.
%
%   SC = CIPIC2SOFAC(CI) allows az, el, and rad to be specified as a 3-
%   column matrix [az,el,rad]. SC is then the 3-column matrix [x,y,z].
%
%   See also SOFAC2CIPIC, SOFAC2SOFAS.

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
    rad = az(:,3);
    el = az(:,2);
    az = az(:,1);
end

az = shiftdim(az);
el = shiftdim(el);
rad = shiftdim(rad);

x = rad.*cosd(-az).*cosd(el);
y = rad.*sind(-az);
z = rad.*cosd(-az).*sind(el);

if singleArgFlag
    x = [x y z];
end

end

