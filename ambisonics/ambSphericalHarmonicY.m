function val = ambSphericalHarmonicY(l, m, r, ambNorm)
%AMBSPHERICALHARMONICY Real-valued spherical harmonic function.
%   Y = AMBSPHERICALHARMONICY(L,M,R,AMBNORM) computes Y, the real-valued
%   spherical harmonic of order L and degree M used in Ambisonics, with
%   normalization convention AMBNORM.
%
%   R may be a P-by-3 matrix of directions, where each row is a Cartesian
%   vector. In this case, Y will be a P-by-1 vector.
%
%   See also AMBNORMALIZATION.

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

%   References:
%     [1] Zotter (2009) Analysis and Synthesis of Sound-Radiation with
%         Spherical Arrays.

% coder.extrinsic('warning')

% Needs at least 3 input arguments
if nargin < 3
    error('Not enough input arguments.');
end

% Uses N3D normalization by default
if nargin < 4
    ambNorm = 'N3D';
end

if (l >= 0) && (abs(m) <= l)
    if isvector(r)
        [AZIM,ELEV,~] = cart2sph(r(1),r(2),r(3));
    else
        [AZIM,ELEV,~] = cart2sph(r(:,1),r(:,2),r(:,3));
    end
    Nlm = ambNormalization(l, abs(m), ambNorm);
    Pl = legendre(l, sin(ELEV));
    Plm = Pl(abs(m) + 1,:).';
%     Plm = legendrePnm(l, abs(m), sin(ELEV));
    
    if m >= 0
        Tm = cos(m * AZIM);
    elseif m < 0
        Tm = sin(abs(m) * AZIM);
    else
        Tm = 0;
    end
    val = Nlm*Plm.*Tm;
else
    warning('Invalid order and degree.');
    val = 0;
end

end