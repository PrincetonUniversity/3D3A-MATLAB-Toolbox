function s = estimateSourceDirection(A_rec, ambNorm, angleBin, force2D)
%ESTIMATESOURCEDIRECTION Source direction estimation via intensity vector.
%   S = ESTIMATESOURCEDIRECTION(A) estimates the source direction vector S
%   given frequency-domain ambisonics signals A.
%
%   S = ESTIMATESOURCEDIRECTION(A,AMBNORM) additionally specifies the
%   ambisonics normalization AMBNORM. By default, N3D is assumed.
%
%   S = ESTIMATESOURCEDIRECTION(A,AMBNORM,BIN) additionally specifies the
%   angular BIN (in degrees) to use for computing the most common
%   direction.
%
%   S = ESTIMATESOURCEDIRECTION(A,AMBNORM,BIN,FORCE2D) optionally forces
%   the estimated source direction to be in the horizontal plane.

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2018 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   ==============================================================================

if nargin < 2 || isempty(ambNorm)
    ambNorm = 'N3D';
end

if nargin < 3 || isempty(angleBin)
    angleBin = 1; % degrees
end

if nargin < 4 || isempty(force2D)
    force2D = false;
end

[rI, ~] = merimaa2005(A_rec, 'ambNorm', ambNorm, 'normalize');
[azim, elev, ~] = cart2sph(rI(:,1),rI(:,2),rI(:,3));

azim(abs(((180/pi)*azim)-360) <= angleBin/2) = 0; % set angles within half a bin of 360 to zero.
theta = mode(angleBin*round((180/pi)*(azim/angleBin)))*pi/180;

if force2D
    phi = 0;
else
    phi = mode(angleBin*round((180/pi)*(elev/angleBin)))*pi/180;
end

[sx,sy,sz] = sph2cart(theta,phi,1);
s = [sx,sy,sz];

end