function psi = getPointSourcePotential(rs,k,r,varargin)
%GETPOINTSOURCEPOTENTIAL Potential due to a point source.
%   PSI = GETPOINTSOURCEPOTENTIAL(RS,K,R) computes the potential field at R
%   due to a point source at RS for an angular wavenumber K. The position 
%   vectors RS and R must be specified in meters in Cartesian coordinates.
%   RS must be a length-3 vector while R must be an M-by-3 matrix where M
%   >= 1. The wavenumber, K, must be specified in rad/m and may be a scalar 
%   or vector of length N. PSI will have dimensions N-by-M is K is a column
%   vector and M-by-N if K is a row vector. The GETPRESSURE function may be 
%   used to compute the pressure from the returned PSI value.
%
%   PSI = GETPOINTSOURCEPOTENTIAL(RS,K,R,TD) optionally allows a scalar 
%   time delay, TD, to be specified in seconds.
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

validateattributes(rs,{'numeric'},{'vector','real','finite','numel',3},...
    'getPointSourcePotential','RS',1);
validateattributes(k,{'numeric'},{'vector','real','finite',...
    'nonnegative'},'getPointSourcePotential','K',2);
validateattributes(r,{'numeric'},{'2d','real','finite','size',[NaN,3]},...
    'getPointSourcePotential','R',3);

% For backwards compatibility
if isrow(k)
    tpFlag = false;
else
    tpFlag = true;
end

rs = shiftdim(rs).'; % Force rs to be a row vector
k = shiftdim(k).'; % Force k to be a row vector
kLen = length(k);
numPos = size(r,1);

if nargin == 4
    c = getSoundSpeed();
    validateattributes(varargin{1},{'numeric'},{'scalar','real',...
        'finite'},'getPointSourcePotential','TD',4);
    phaseDelay = repmat(exp(1i*k*c*varargin{1}),numPos,1);
else
    phaseDelay = 1; % No delay by default
end

dr = r-repmat(rs,numPos,1);
d = computeVectorNorm(dr,2,2);
preTerm = repmat(1./d,1,kLen);

% For backwards compatibility
if tpFlag
    psi = (preTerm.*exp(1i*d*k).*phaseDelay).';
else
    psi = preTerm.*exp(1i*d*k).*phaseDelay;
end

end
