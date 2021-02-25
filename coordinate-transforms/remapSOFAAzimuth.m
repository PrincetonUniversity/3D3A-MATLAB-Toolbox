function azOut = remapSOFAAzimuth(azIn)
%REMAPSOFAAZIMUTH Modify range of azimuth specified in SOFA spherical
%coordinates.
%   B = REMAPSOFAAZIMUTH(A) takes azimuth, A, in degrees, and returns the
%   remapped azimuth, B, in degrees. A may be a scalar, vector, or matrix. 
%   B will have the same dimensions as A.
%       If A is in the range [0,360), B will be in the range (-180,180].
%
%       If A is in the range (-180,180], B will be in the range [0,360).
%
%   In all other cases, an error message will be triggered.

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

narginchk(1,1)

validateattributes(azIn,{'numeric'},{'2d','nonempty','finite','real',...
    '>',-180,'<',360},'remapSOFAAzimuth','A',1)

if any(any(azIn > 180)) && ~any(any(azIn < 0))
    azOut = azIn; % Initialize azOut
    remapIndxs = azIn > 180; % Find indices where azimuth exceeds 180
    azOut(remapIndxs) = azIn(remapIndxs)-360; % Remap azimuth
elseif any(any(azIn < 0)) && ~any(any(azIn > 180))
    azOut = mod(mod(azIn,360),360); % Remap azimuth
elseif any(any(azIn < 0)) && any(any(azIn > 180))
    error('Invalid input. A must be in the range (-180,180] or [0,360).')
else
    azOut = azIn;
end

end
