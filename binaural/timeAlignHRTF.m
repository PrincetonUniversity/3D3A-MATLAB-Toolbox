function HOut = timeAlignHRTF(HIn,SP,F,varargin)
%TIMEALIGNHRTF Time align HRTF using a rigid-sphere model.
%   HOut = TIMEALIGNHRTF(HIn,SP,F) attempts to modify the linear phase
%   component of an input HRTF, HIn, based on a rigid-sphere model (i.e., 
%   assuming the head is a rigid-sphere) so that the HRTFs corresponding to
%   different source directions all end up having approximately the same 
%   linear phase component (i.e., are aligned in time). HIn must be an 
%   M-by-N matrix of frequency responses representing the HRTF for a single 
%   ear and for source positions, SP, which must be an N-by-3 matrix 
%   specifying positions in SOFA cartesian coordinates. F must be a length 
%   M vector of frequencies, in Hz, at which the HRTF is specified. It is 
%   assumed that the radius of the rigid-sphere is A = 8.75 cm, and that 
%   the location of the ear on the sphere is [90,0,A] as specified in SOFA 
%   spherical coordinates. This corresponds to the left ear (i.e., it is 
%   assumed that the specified HIn corresponds to a left ear HRTF). These 
%   values can be changed by specifying optional parameters (see below). 
%   HOut will have the same dimensions as HIn.
%
%   HOut = TIMEALIGNHRTF(HIn,SP,F,A) optionally allows the radius, A, of 
%   the rigid sphere to be specified in meters. If A is not specified, or 
%   is specified as [], a default value of 0.0875 m is used. A must be a 
%   scalar.
%
%   HOut = TIMEALIGNHRTF(HIn,SP,F,A,EP) optionally allows the direction, 
%   EP, of the ear to be specified. EP must be specified as [az,el], where 
%   az and el are azimuth and elevation in degrees following the SOFA 
%   spherical coordinate system conventions. az can vary from 0 to 360 and 
%   el from -90 to 90. If EP is not specified or is specifed as [], a 
%   default value of EP = [90,0] is used, and corresponds to the left ear 
%   located on the horizontal plane. EP must be a 1-by-2 vector.
%
%   HOut = TIMEALIGNHRTF(HIn,SP,F,A,EP,R) optionally allows source 
%   distance, in meters, to be explicitly specified. Specifying R overrides 
%   source distance computed directly from SP. If R is not specified or is 
%   specified as [], the default value is computed from SP. R can be a 
%   scalar, including inf (corresponding to the theoretical far-field), or 
%   an N-by-1 vector.
%
%   HOut = TIMEALIGNHRTF(HIn,SP,F,A,EP,R,METHOD) optionally allows 
%   specification of the method of computing the modification to be applied 
%   to the linear phase component. METHOD can take one of the following 
%   values:
%       1. 'Ben-HurEtAl2019' - use Eq. (13) provided by Ben-Hur et al. [1].
%
%       2. 'Ben-HurEtAl2019nf' - extend Eq. (13) provided by Ben-Hur 
%       et al. [1]. to also apply to near-field source distances.
%
%       3. 'SridharChoueiri2021a' - use Eqs. (7) and (8) provided by 
%       Sridhar and Choueiri [2].
%
%       4. 'SridharChoueiri2021b' - use Eqs. (9) and (10) provided by 
%       Sridhar and Choueiri [2]. This is the default.

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
%   Copyright (c) 2021 Princeton University
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

%   Refs:
%       [1]. Ben-Hur, Z., Alon, D. L., Mehra, R., and Rafaely, B.,
%       "Efficient Representation and Sparse Sampling of Head-Related 
%       Transfer Functions Using Phase- Correction Based on Ear Alignment,"
%       IEEE/ACM Transactions on Audio, Speech, and Language Pro- cessing, 
%       27(12), pp. 2249?2262, 2019.
%
%       [2]. Sridhar, R. and Choueiri, E., "Validation of the Minimum-Phase 
%       Nature of the Rigid-Sphere Head-Related Transfer Function with 
%       Applications in Onset Estimation," J. Audio Eng. Soc, 2021, Under 
%       Review.

narginchk(3,7);

numOpts = length(varargin);
switch numOpts
    case 1 % Only A is specified
        A = varargin{1};        
        EP = [90,0];
        R = [];
        computeRFlag = true;
        METHOD = 'SridharChoueiri2021b';
    case 2 % Only A and EP are specified
        A = varargin{1};
        EP = varargin{2};
        R = [];
        computeRFlag = true;
        METHOD = 'SridharChoueiri2021b';
    case 3 % Only A, EP and R are specified
        A = varargin{1};
        EP = varargin{2};
        R = varargin{3};
        computeRFlag = false;
        METHOD = 'SridharChoueiri2021b';
    case 4 % All are specified
        A = varargin{1};
        EP = varargin{2};
        R = varargin{3};
        METHOD = varargin{4};
        computeRFlag = false;
    otherwise % None are specified
        A = 0.0875;
        EP = [90,0];
        R = [];
        computeRFlag = true;
        METHOD = 'SridharChoueiri2021b';
end

% Apply default values if necessary
if isempty(A)
    A = 0.0875;
end
if isempty(EP)
    EP = [90,0];
end
if isempty(R)
    computeRFlag = true;
end

% Compute angles of incidence, theta
numPos = size(SP,1);
SD = sofaC2sofaS(SP);
[xS,yS,zS] = sofaS2sofaC(SD(:,1),SD(:,2),ones(numPos,1));
[xE,yE,zE] = sofaS2sofaC(EP(1),EP(2),1);
theta = getCentralAngle([xS,yS,zS],repmat([xE,yE,zE],numPos,1));

% Compute non-dimensional quantities and other required variables
c = getSoundSpeed();
F = shiftdim(F); % Force F to be a column vector
mu = 2*pi*F*A/c;
if computeRFlag
    R = SD(:,3);
else
    if isscalar(R)
        R = R*ones(numPos,1);
    end
end
rho = R/A;
if any(rho <= 1)
    error(['Specified SP incompatible with A = %f m. A must be smaller',...
        ' than the source distance.'],A)
end
x = cosd(theta);
g = sqrt(rho.^2-(2*rho.*x)+1);

% Compute and remove linear phase component
HOut = zeros(size(HIn)); % Initialization only
switch lower(METHOD)
    case 'ben-huretal2019'
        HOut = HIn.*exp(-1i*mu*(x.'));
    case 'ben-huretal2019nf'
        for ii = 1:numPos
            if rho(ii) == inf
                HOut(:,ii) = HIn(:,ii).*exp(-1i*mu*x(ii));
            else
                HOut(:,ii) = HIn(:,ii).*exp(1i*mu*g(ii)).*exp(-1i*mu*...
                    rho(ii));
            end
        end
    case'sridharchoueiri2021a'
        for ii = 1:numPos
            if rho(ii) == inf
                if theta(ii) <= 90
                    delay = mu*x(ii);
                else
                    delay = -mu*deg2rad(theta(ii)-90);
                end
            else
                theta0 = acosd(1/rho(ii));
                hatG = sqrt(rho(ii)^2-1);
                if theta(ii) < theta0
                    delay = -mu*g(ii) + mu*rho(ii);
                else
                    delay = -mu*(hatG + deg2rad(theta(ii)-theta0)) + ...
                        mu*rho(ii);
                end
            end
            HOut(:,ii) = HIn(:,ii).*exp(-1i*delay);
        end
    case 'sridharchoueiri2021b'
        mu_cw = 2*pi*max(F)*A/c;
        corFac = 1/(1+0.5094*(2*mu_cw^2)^(-1/3));
        
        for ii = 1:numPos
            if rho(ii) == inf
                if theta(ii) <= 90
                    delay = mu*x(ii);
                else
                    delay = -mu*deg2rad(theta(ii)-90)/corFac;
                end
            else
                theta0 = acosd(1/rho(ii));
                hatG = sqrt(rho(ii)^2-1);
                if theta(ii) < theta0
                    delay = -mu*g(ii) + mu*rho(ii);
                else
                    delay = -(mu*hatG + mu*deg2rad(theta(ii)-theta0)/...
                        corFac) + mu*rho(ii);
                end
            end
            HOut(:,ii) = HIn(:,ii).*exp(-1i*delay);
        end
    otherwise
        error('Invalid METHOD specification.')
end

end
