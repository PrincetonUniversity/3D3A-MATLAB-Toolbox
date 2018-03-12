function T = AmbiNav_translation(Li, Lo, d, kVec)
%AMBINAV_TRANSLATION Ambisonics translation coefficients matrix.
%   T = AMBINAV_TRANSLATION(LI,LO,D,K) computes the ambisonic translation
%   coefficients matrix T, for input ambisonics order LI, output order LO,
%   translation position vector D (given in Cartesian coordinates), and for
%   angular wavenumber K. K may be a vector, in which case T is (LO+1)^2-by
%   -(LI+1)^2-by-LENGTH(K). The N3D ambisonics normalization convention is
%   assumed.

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
%   Copyright (c) 2017 Princeton University
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
%     [1] Gumerov and Duraiswami (2005) Fast Multipole Methods for the
%         Helmholtz Equation in Three Dimensions.
%     [2] Zotter (2009) Analysis and Synthesis of Sound-Radiation with
%         Spherical Arrays.

narginchk(4,4);

if numel(d) == 3
    [AZIM,ELEV,R] = cart2sph(d(1),d(2),d(3));
else
    error('Translation vector D should have three elements.');
end

kLen = length(kVec);

Ni = (Li + 1)^2;
No = (Lo + 1)^2;

if R == 0
    QzL = eye(No);
    QzR = eye(Ni);
else
    QzL = AmbiNav_zRotation(AZIM, ELEV, Lo);
    QzR = AmbiNav_zRotation(AZIM, ELEV, Li);
end

T = zeros(No,Ni,kLen);
Tz = AmbiNav_zTranslation(kVec*R, max([Li, Lo]));
for kk = 1:kLen
    if kVec(kk)*R < f2k(10)*0.001
        T(:,:,kk) = eye(No,Ni);
    else
        T(:,:,kk) = QzL*Tz(1:No,1:Ni,kk)/QzR;
    end
end

end