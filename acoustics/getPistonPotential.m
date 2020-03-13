function psi = getPistonPotential(r0,k,a,r,varargin)
%getPistonPotential Potential due to a circular piston.
%   psi = getPistonPotential(r0,k,a,r) computes the potential field at r
%   due to a baffled piston source of radius, a, at r0 for an angular 
%   wavenumber k. The position vectors r0 and r should be specified in 
%   meters in SOFA Cartesian coordinates and must each be P-by-3 matrices, 
%   where P corresponds to the number of positions. The piston radius, a, 
%   should be a scalar specified in meters, while the wavenumber, k, may be 
%   a scalar or vector of length N in units or rad/m. psi will dimensions 
%   of N-by-P. The GETPRESSURE function may be used to compute the pressure 
%   from the returned psi value.
%
%   psi = getPistonPotential(r0,k,a,r,u0n) optionally allows a scalar 
%   piston acceleration, u0n, to be specified in m/s^2. If u0n is not
%   specified, or is specified as [], a default value of 1/(1.21*a^2) is 
%   used.
%
%   psi = getPistonPotential(r0,k,a,r,u0n,t0) optionally allows a scalar
%   time delay, t0, to be specified in seconds.
%
%   See also GETPOINTSOURCEPOTENTIAL, GETPRESSURE.

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

% Check number of inputs
narginchk(4,6);

k = shiftdim(k).';
kLen = length(k);
numPos = size(r,1);
if size(r0,1) ~= size(r,1)
    r0 = repmat(r0,numPos,1);
end

if nargin < 6
    phaseDelay = 1; % No delay by default
else
    c = getSoundSpeed();
    phaseDelay = repmat(exp(1i*k*c*varargin{2}),numPos,1);
end

if nargin < 5
    u0n = 1/(1.21*a^2);
else
    if isempty(varargin{1})
        u0n = 1/(1.21*a^2);
    else
        u0n = varargin{1};
    end
end

dr = r-r0;
d = vecnorm(dr,2,2);
st = dr(:,2)./d;
preTerm = repmat(1.21*u0n*a^2./d,1,kLen);

psi = (preTerm.*exp(1i*d*k).*(besselj(1,st*k*a)./(st*k*a)).*phaseDelay).';

end
