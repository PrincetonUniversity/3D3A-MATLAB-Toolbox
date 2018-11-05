function [az,el,rad] = sofaC2cipicI(x,y,z,FLAG)
%SOFAC2CIPICI Convert SOFA cartesian coordinates to CIPIC interaural 
%coordinates.
%   [az,el,rad] = SOFAC2CIPICI(x,y,z) converts from SOFA cartesian  
%   coordinates to CIPIC interaural coordinates. az and el are specified in
%   degrees. x, y, and z may be specified as scalars or vectors. If 
%   vectors, they must have the same length. az, el, and rad are either 
%   scalars or column vectors.
%
%   CI = SOFAC2CIPICI(SC) allows x, y, and z to be specified as a 3-column 
%   matrix [x,y,z]. CI is then the 3-column matrix [az,el,rad].
%
%   __ = SOFAC2CIPICI(__,'flipAz') optionally allows the sign convention of
%   azimuth to be flipped such that positions corresponding to positive
%   values of 'y' in SOFA cartesian coordinates correspond to positive
%   azimuths (unlike in the CIPIC coordinate system, where these positions
%   correspond to negative azimuths).
%
%   See also CIPICI2SOFAC, SOFAS2SOFAC.

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

narginchk(1,4);

if nargin <= 2 && nargout <= 1
    singleArgFlag = true;
    if nargin == 2
        FLAG = y;
    else
        FLAG = 'noFlipAz';
    end
    z = x(:,3);
    y = x(:,2);
    x = x(:,1);
else
    singleArgFlag = false;
    if nargin < 4
        FLAG = 'noFlipAz';
    end
    z = shiftdim(z);
    y = shiftdim(y);
    x = shiftdim(x);
end

rad = sqrt(x.^2 + y.^2 + z.^2);
I = rad ~= 0;
rad = rad(I);
az = zeros(size(rad));

if isscalar(x)
    x = ones(size(rad))*x;
end
if isscalar(y)
    y = ones(size(rad))*y;
end
if isscalar(z)
    z = ones(size(rad))*z;
end

switch lower(FLAG)
    case {'flipaz'}
        az(I) = -asind(-y(I)./rad);
    case {'noflipaz'}
        az(I) = asind(-y(I)./rad);
    otherwise
        error('Invalid FLAG specification')
end
el = mod(atan2d(z(I),x(I)),360);
el(el >= 270) = el(el >= 270)-360;

if singleArgFlag
    az = [az el rad];
end

end
