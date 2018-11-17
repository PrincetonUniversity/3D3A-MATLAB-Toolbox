function val = ambHelmholtzSolutionR(l,m,k,r,ambNorm)
%AMBHELMHOLTZSOLUTIONR Regular solution to 3D Helmholtz equation.
%   V = AMBHELMHOLTZSOLUTIONR(L,M,K,R) computes V, the regular solution to
%   the 3D Helmholtz equation given by:
%
%       V = 4 * PI * (-I)^L * j_L(K*NORM(R)) * Y_L^M(R) / ||Y_L^M||^2
%
%   where K is the angular wavenumber, R is the position (given in
%   Cartesian coordinates), j_L is the spherical Bessel function of order
%   L, Y_L^M is the real-valued spherical harmonic of degree L and order M,
%   and ||Y_L^M||^2 is its squared-norm.
%
%   K may be a vector and R may be a P-by-3 matrix, in which case V will be
%   a LENGTH(K)-by-P matrix.
%
%   V = AMBHELMHOLTZSOLUTIONR(L,M,K,R,AMBNORM) additionally specifies the
%   ambisonics normalization convention to use. By default, N3D is assumed.
%
%   See also SPHERICALBESSELJ, AMBSPHERICALHARMONICY, AMBNORMSQUARED.

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
%     [1] Gumerov and Duraiswami (2005) Fast Multipole Methods for the
%         Helmholtz Equation in Three Dimensions.
%     [2] Zotter (2009) Analysis and Synthesis of Sound-Radiation with
%         Spherical Arrays.

narginchk(4,5);

% Uses N3D normalization by default
if nargin < 5 || isempty(ambNorm)
    ambNorm = 'N3D';
end

if (l >= 0) && (abs(m) <= l)
    coeff = 4 * pi * ((-1i)^l);
    jl = sphericalBesselJ(l,k(:)*(sqrt(dot(r,r,2)).'));
    Ylm = ambSphericalHarmonicY(l,m,r,ambNorm);
    val = coeff * jl * diag(Ylm) / ambNormSquared(l,ambNorm);
else
    error('Invalid order and degree.');
end

end