function Nlm = ambNormalization(l, m, ambNorm)
%AMBNORMALIZATION Ambisonics spherical harmonic normalization factor.
%   N = AMBNORMALIZATION(L,M,AMBNORM) returns the normalization factor
%   for the real-valued spherical harmonics used in Ambisonics.
%
%   See also AMBSPHERICALHARMONICY.

%   =======================================================================
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

%   References:
%     [1] Zotter (2009) Analysis and Synthesis of Sound-Radiation with
%         Spherical Arrays.

% Needs at least 2 input arguments
narginchk(2,3);

% Uses N3D normalization by default
if nargin < 3
    ambNorm = 'N3D';
end

switch lower(ambNorm)
    case 'none' % For compatibility with ambisonics.ch website
        Nlm = ((-1)^m)*sqrt((2*l+1).*(2-(~m))).*sqrt(factorial(l-m)./factorial(l+m));
    case 'sn3d' % Not confirmed
        Nlm = ((-1)^m)*sqrt((2-(~m))/(4*pi)).*sqrt(factorial(l-m)./factorial(l+m));
    case 'n3d' % Eq. (31)
        Nlm = ((-1)^m)*sqrt((2*l+1).*(2-(~m))/(4*pi)).*sqrt(factorial(l-m)./factorial(l+m));
    otherwise
        error('Unknown normalization convention.');
end

end