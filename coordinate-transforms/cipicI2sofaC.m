function [x,y,z] = cipicI2sofaC(az,el,rad,FLAG)
%CIPICI2SOFAC Convert CIPIC interaural coordinates to SOFA cartesian
%coordinates.
%   [x,y,z] = CIPICI2SOFAC(az,el,rad) converts from CIPIC interaural  
%   coordinates to SOFA cartesian coordinates. Input azimuth (az) and 
%   elevation (el) must be specified in degrees. az, el, and rad may be 
%   specified as scalars or vectors. If vectors, they must have the same 
%   length. x, y, and z have the same dimensions as az, el, and rad.
%
%   SC = CIPICI2SOFAC(CI) allows az, el, and rad to be specified as a 3-
%   column matrix [az,el,rad]. SC is then the 3-column matrix [x,y,z].
%
%   ___ = CIPICI2SOFAC(...,'flipAz') optionally indicates that the sign 
%   convention of azimuth is flipped such that positions corresponding to 
%   positive values of 'y' in SOFA cartesian coordinates correspond to 
%   positive azimuths (unlike in the CIPIC coordinate system, where these 
%   positions correspond to negative azimuths).
%
%   See also SOFAC2CIPICI, SOFAC2SOFAS.

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
    
    if size(az,2) ~= 3
        error('Invalid input dimensions')
    end
    
    if nargin == 2
        FLAG = el;
        if ~ischar(FLAG)
            error('FLAG must be of type char')
        end
        
        if strcmpi(FLAG,'flipAz')
            rad = az(:,3);
            el = az(:,2);
            az = -az(:,1);
        else
            error('Invalid FLAG specification')
        end
    else
        rad = az(:,3);
        el = az(:,2);
        az = az(:,1);
    end
else
    singleArgFlag = false;
    
    if ~isvector(az) || ~isvector(el) || ~isvector(rad)
        error(['az, el, and rad must be scalars or vectors when ',...
            'specified independently'])
    end
    
    if (length(az) ~= length(el)) || (length(el) ~= length(rad)) || ...
            (length(rad) ~= length(az))
        error(['az, el, and rad must have the same length when ',...
            'specified independently'])
    end
    
    if nargin < 4
        rad = shiftdim(rad);
        el = shiftdim(el);
        az = shiftdim(az);
    else
        if ~ischar(FLAG)
            error('FLAG must be of type char')
        end
        
        if strcmpi(FLAG,'flipAz')
            rad = shiftdim(rad);
            el = shiftdim(el);
            az = shiftdim(-az);
        else
            error('Invalid FLAG specification')
        end
    end
end

if (max(az) > 90) || (min(az) < -90)
    error('Azimuths must range from -90 to 90')
end

if (max(el) > 270) || (min(el) < -90)
    error('Elevations must range from -90 to 270')
end

x = rad.*cosd(-az).*cosd(el);
y = rad.*sind(-az);
z = rad.*cosd(-az).*sind(el);

if singleArgFlag
    x = [x y z];
end

end
