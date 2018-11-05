function [az,el,rad] = sofaC2sofaS(x,y,z)
%SOFAC2SOFAS Convert SOFA cartesian coordinates to SOFA spherical 
%coordinates.
%   [az,el,rad] = SOFAC2SOFAS(x,y,z) converts from SOFA cartesian 
%   coordinates to SOFA spherical coordinates. x, y, and z may be scalars 
%   or vectors. az and el are specified in degrees. az, el, and rad are 
%   either scalars, or column vectors with the same length as x, y, and z.
%
%   SS = SOFAC2SOFAS(SC) allows x, y, and z to be specified as a 3-column
%   matrix [x,y,z]. SS is then the 3-column matrix [az,el,rad].
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
    
    validateattributes(x,{'double'},{'2d','nonempty','nonnan','finite',...
        'real','size',[NaN,3]},'sofaC2sofaS','SC',1)
    
    z = x(:,3);
    y = x(:,2);
    x = x(:,1);
elseif nargin == 3
    
    validateattributes(x,{'double'},{'vector','nonempty','nonnan',...
        'finite','real'},'sofaC2sofaS','x',1)
    validateattributes(y,{'double'},{'vector','nonempty','nonnan',...
        'finite','real','numel',numel(x)},'sofaC2sofaS','y',2)
    validateattributes(z,{'double'},{'vector','nonempty','nonnan',...
        'finite','real','numel',numel(x)},'sofaC2sofaS','z',3)
    
    x = shiftdim(x);
    y = shiftdim(y);
    z = shiftdim(z);
else
    error('2 inputs provided when either 1 or 3 inputs are required.')
end

rad = sqrt(x.^2 + y.^2 + z.^2);

if rad == 0
    error('x, y, and z cannot all be zero.')
end

az = mod(atan2d(y,x),360);
el = asind(z./rad);

if singleArgFlag
    az = [az el rad];
end

end
