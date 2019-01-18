function val = sphericalBesselJ(n,x)
%SPHERICALBESSELJ Spherical Bessel function.
%   J = SPHERICALBESSELJ(N,X) returns the value J of the spherical Bessel
%   function of the first kind for order N and with argument X.
%
%   N and X may be matrices, in which case they must have the same size or
%   else one must be a scalar.
%
%   See also BESSELJ.

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
%   =======================================================================

if ~(all(size(x)==size(n)) || isscalar(x) || isscalar(n))
    error('Inputs N and X must be the same size or else one must be a scalar.');
end

coeff = sqrt(pi./(2*x));
sgn = 2*(x>=0)-1; % Signum-like function; equal to +1 at 0
val = coeff.*sgn.*besselj(n+0.5,x);
if isscalar(n)
    val(x==0) = +~n; % gives 1 if n==0; gives 0 if n~=0
elseif isscalar(x) && (x==0)
    val = +~n;
else
    val(x==0) = +~n(x==0);
end

end