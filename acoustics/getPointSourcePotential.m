function psi = getPointSourcePotential(r0,k,r,varargin)
%getPointSourcePotential Potential due to a point source.
%   psi = getPointSourcePotential(r0,k,r) computes the potential field at r
%   due to a point source at r0 for an angular wavenumber k. The position 
%   vectors r0 and r should be specified in meters in SOFA Cartesian 
%   coordinates and must each be three-element vectors. The wavenumber, k, 
%   may be a scalar or vector in units or rad/m. psi will have the same 
%   dimensions as k. The GETPRESSURE function may be used to compute the
%   pressure from the returned psi value.
%
%   psi = getPointSourcePotential(r0,k,r,t0) optionally allows a scalar 
%   time delay, t0, to be specified in seconds.
%
%   See also GETPISTONPOTENTIAL, GETPRESSURE.

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
narginchk(3,4);

if nargin == 4
    c = getSoundSpeed();
    phaseDelay = exp(1i*k*c*varargin{1});
else
    phaseDelay = 1; % No delay by default
end

d = norm(r - r0);
psi = (exp(1i*k*d)/d).*phaseDelay;

end
