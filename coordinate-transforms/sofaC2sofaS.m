function [az,el,rad] = sofaC2sofaS(x,y,z)
%SOFAC2SOFAS Convert SOFA cartesian coordinates to SOFA spherical 
%coordinates.
%   [az,el,rad] = SOFAC2SOFAS(x,y,z) converts from SOFA cartesian to
%   SOFA spherical coordinates. x, y, and z may be scalars or vectors. az
%   and el are specified in degrees. az, el, and rad may be scalars or
%   column vectors.
%
%   [RS] = SOFAC2SOFAS(RC) allows x, y, and z to be specified as a 3-column
%   matrix [x,y,z]. RS is then the 3-column matrix [az,el,rad].
%
%   See also SOFAS2SOFAC.

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
    z = x(:,3);
    y = x(:,2);
    x = x(:,1);
end

x = shiftdim(x);
y = shiftdim(y);
z = shiftdim(z);

rad = sqrt(x.^2 + y.^2 + z.^2);
I = rad ~= 0;
rad = rad(I);

az = mod(atan2d(y(I),x(I)),360);
el = asind(z(I)./rad);

if singleArgFlag
    az = [az el rad];
end

end
