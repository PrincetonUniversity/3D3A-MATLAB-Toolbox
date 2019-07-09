function [H,N,thVec] = computeSphereHRTF(a,r,theta,f,METHOD)
%COMPUTESPHEREHRTF Analytically computed HRTF for a sphere.
%   [H,N,T] = COMPUTESPHEREHRTF(A,R,THETA,F) analytically computes the HRTF
%   value, H, for a sphere of radius A (in meters), source at distance R
%   (in meters) with angle of incidence THETA (in degrees), and for
%   frequency F (in Hz) using the method described by Sridhar and Choueiri
%   [1] with a numerical precision value of inf (i.e., maximum precision).
%   All inputs must be real-valued scalars. THETA can take values in the 
%   range [0,360). Also returned are the computed order, N, used to compute 
%   H, and a 2-by-1 vector, T, of corresponding threshold values for use 
%   with the HRTF computation methods described by Cooper and Bauck [2], 
%   and Duda and Martens [3], in that order. The GETCENTRALANGLE function 
%   may be useful for computing THETA. 
%
%   [H,N,T] = COMPUTESPHEREHRTF(A,R,THETA,F,METHOD) optionally specifies
%   the METHOD to use to compute the HRTF. METHOD must be a 2-element cell 
%   array and can take the following values:
%       1. {'sridharchoueiri2019',P} where P can be any finite integer or 
%       inf. The value of precision, P, refers to the number of decimal 
%       places to use when determining the order, N, of the calculation, 
%       with P = inf corresponding to maximum precision. If METHOD is not 
%       specified, {'sridharchoueiri2019',inf} is used.
%       2. {'fixedn',N} where N is the order and must be a non-negative
%       integer. In this case, the second output, N, will be the same as 
%       this specified input order.
%       3. {'cooperbauck1980',TH} where TH can be any finite, real number,
%       and refers to the threshold for determining the order, N, of the
%       calculation. In the Fortran code published by Cooper and Bauck [2], 
%       TH is set to 0.001. In this case, the specified value of TH will be
%       also be returned as the first element of the output vector, T.
%       4. {'dudamartens1998',TH} where TH is as defined in 3, above. In
%       this case, the specified value of TH will be also be returned as 
%       the second element of the output vector, T.
%       5. {'exact',P} where P is as defined in 1, above.
%
%   See also GETCENTRALANGLE, COMPUTESPHEREITD.

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
%   Copyright (c) 2019 Princeton University
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

% Refs:
%   [1]. Sridhar and Choueiri (to be published) - Computational Analysis 
%   and Minimum-Phase Characteristics of the Head-Related Transfer 
%   Functions of a Rigid Sphere.
%   [2]. Cooper and Bauck (1980) - On Acoustical Specification of Natural 
%   Stereo Imaging.
%   [3]. Duda and Martens (1998) - Range dependence of the response of a 
%   spherical head model.

narginchk(4,5);

if nargin < 5
    METHOD{1} = 'sridharchoueiri2019';
    METHOD{2} = inf; % precision/threshold, depending on METHOD{1}
                     % Note: threshold = 10^(-precision)
end

% Check inputs
validateattributes(a,{'double'},{'scalar','nonempty','nonnan','finite',...
    'positive','real'},'computeSphereHRTF','A',1);
validateattributes(r,{'double'},{'scalar','nonempty','nonnan'},...
    'computeSphereHRTF','R',2);
validateattributes(theta,{'double'},{'scalar','nonempty','nonnan',...
    'finite','real','nonnegative','<',360},'computeSphereHRTF','THETA',3);
validateattributes(f,{'double'},{'scalar','nonempty','nonnan','finite',...
    'nonnegative','real'},'computeSphereHRTF','F',4);
validateattributes(METHOD{1},{'char'},{'scalartext','nonempty'},...
    'computeSphereHRTF','first element of METHOD',5);
validateattributes(METHOD{2},{'double'},{'scalar','nonempty','nonnan'},...
    'computeSphereHRTF','second element of METHOD',5);

% Initialize common variables
c = getSoundSpeed();
mu = (2*pi*f*a)/c;
thVec = zeros(2,1);
x = cosd(theta);
rho = r/a;

% Main calculation
if mu == 0
    H = 1;
    N = 0;
    thVec = ones(2,1);
else
    % Perform preliminary calculations
    % Note: hf is the variable storing the value of the spherical-Hankel 
    % function (of the first kind) with mu as the argument. P is the 
    % variable storing the value of the Legendre polynomial.
    
    hf = zeros(3,1); % Initialize hf to store 3 most recent terms
    P = zeros(3,1); % Initialize P to store 3 most recent terms
    hf(1) = exp(1i*mu)/(1i*mu); % m = 0
    hf(2) = hf(1)*((1/mu)-1i); % m = 1
    hf(3) = computeHX(hf(1:2),2,mu); % m = 2
    P(1) = 1; % m = 0
    P(2) = x; % m = 1
    
    % Initialize sums
    psiPS_old = 0;
    psiPS = 0;
    switch lower(METHOD{1})
        case 'fixedn'
            N = METHOD{2}; % extract order from input
            if rho == inf
                m = 0;
                dh = computeDH(hf(1:2),m,mu); % m = 0
                Am = computeAm(m,dh); % m = 0
                psiP = conj(Am)*P(1);
                psiPS = psiPS + psiP;
                m = m + 1;
                % Store most recent values
                psiP_old = psiP;
                psiPS_older = psiPS_old;
                psiPS_old = psiPS;
                % Compute new values
                dh = computeDH(hf(2:3),m,mu); % m = 1
                Am = computeAm(m,dh); % m = 1
                psiP = conj(Am)*P(2);
                psiPS = psiPS + psiP;
                for m = 2:N
                    % Store old values
                    psiP_old = psiP;
                    psiPS_older = psiPS_old;
                    psiPS_old = psiPS;
                    % Compute new values
                    hf(1) = hf(2);
                    hf(2) = hf(3);
                    hf(3) = computeHX(hf(1:2),m+1,mu); % m+1 is ok, not m
                    dh = computeDH(hf(2:3),m,mu);
                    Am = computeAm(m,dh);
                    P(3) = computeP(P(1:2),m,x);
                    psiP = conj(Am)*P(3);
                    psiPS = psiPS + psiP;
                    P(1) = P(2);
                    P(2) = P(3);
                end
                H = (1/mu^2)*psiPS;
            else
                mr = mu*rho;
                hn = zeros(3,1); % Init. hn to store 3 most recent terms
                % Note: hn is the variable storing the value of the 
                % spherical-Hankel function (of the first kind) with mr as 
                % the argument.
                
                m = 0;
                dh = computeDH(hf(1:2),m,mu); % m = 0
                hn(1) = exp(1i*mr)/(1i*mr); % m = 0
                Bm = computeBm(m,dh,hn(1)); % m = 0
                psiP = conj(Bm)*P(1);
                psiPS = psiPS + psiP;
                m = m + 1;
                % Store old values
                psiP_old = psiP;
                psiPS_older = psiPS_old;
                psiPS_old = psiPS;
                % Compute new values
                dh = computeDH(hf(2:3),m,mu); % m = 1
                hn(2) = hn(1)*((1/mr)-1i); % m = 1
                Bm = computeBm(m,dh,hn(2)); % m = 1
                psiP = conj(Bm)*P(2);
                psiPS = psiPS + psiP;
                for m = 2:N
                    % Store old values
                    psiP_old = psiP;
                    psiPS_older = psiPS_old;
                    psiPS_old = psiPS;
                    % Compute new values
                    hf(1) = hf(2);
                    hf(2) = hf(3);
                    hf(3) = computeHX(hf(1:2),m+1,mu); % m+1 is ok, not m
                    dh = computeDH(hf(2:3),m,mu);
                    hn(3) = computeHX(hn(1:2),m,mr); % m is ok, not m+1
                    Bm = computeBm(m,dh,hn(3));
                    P(3) = computeP(P(1:2),m,x);
                    psiP = conj(Bm)*P(3);
                    psiPS = psiPS + psiP;
                    P(1) = P(2);
                    P(2) = P(3);
                    hn(1) = hn(2);
                    hn(2) = hn(3);
                end
                H = -(rho/mu)*exp(-1i*mu*rho)*psiPS;
            end
            % Compute equivalent thresholds
            for ii = 1:2
                [~,thVec(ii),~] = evaluateConvergence(psiPS_older,...
                    psiP_old,psiPS_old,psiP,psiPS,0,ii); % 0 - dummy value
            end
        case {'cooperbauck1980','dudamartens1998','sridharchoueiri2019',...
                'exact'}
            switch lower(METHOD{1})
                case 'cooperbauck1980'
                    methodFlag = 1;
                case 'dudamartens1998'
                    methodFlag = 2;
                case 'sridharchoueiri2019'
                    methodFlag = 3;
                case 'exact'
                    methodFlag = 4;
                otherwise
                    error('Invalid METHOD specification.')
            end
            thVal = METHOD{2}; % extract threshold from input
            if rho == inf
                m = 0;
                dh = computeDH(hf(1:2),m,mu); % m = 0
                Am = computeAm(m,dh); % m = 0
                psiP = conj(Am)*P(1);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                m = m + 1;
                % Store old values
                psiP_old = psiP;
                psiPS_older = psiPS_old;
                psiPS_old = psiPS;
                % Compute new values
                dh = computeDH(hf(2:3),m,mu); % m = 1
                Am = computeAm(m,dh); % m = 1
                psiP = conj(Am)*P(2);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                [~,~,lF] = evaluateConvergence(psiPS_older,psiP_old,...
                    psiPS_old,psiP,psiPS,thVal,methodFlag);
                while lF
                    m = m + 1;
                    % Store old values
                    psiP_old = psiP;
                    psiPS_older = psiPS_old;
                    psiPS_old = psiPS;
                    % Compute new values
                    hf(1) = hf(2);
                    hf(2) = hf(3);
                    hf(3) = computeHX(hf(1:2),m+1,mu); % m+1 is ok, not m
                    dh = computeDH(hf(2:3),m,mu);
                    Am = computeAm(m,dh);
                    P(3) = computeP(P(1:2),m,x);
                    psiP = conj(Am)*P(3);
                    psiPS = psiPS + psiP;
                    psiPCS{m+1} = psiPS;
                    P(1) = P(2);
                    P(2) = P(3);
                    [~,~,lF] = evaluateConvergence(psiPS_older,...
                        psiP_old,psiPS_old,psiP,psiPS,thVal,methodFlag);
                end
                H = (1/mu^2)*psiPS;
            else
                mr = mu*rho;
                hn = zeros(3,1); % Init. hn to store 3 most recent terms
                % Note: hn is the variable storing the value of the 
                % spherical-Hankel function (of the first kind) with mr as 
                % the argument.
                
                m = 0;
                dh = computeDH(hf(1:2),m,mu); % m = 0
                hn(1) = exp(1i*mr)/(1i*mr); % m = 0
                Bm = computeBm(m,dh,hn(1)); % m = 0
                psiP = conj(Bm)*P(1);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                m = m + 1;
                % Store old values
                psiP_old = psiP;
                psiPS_older = psiPS_old;
                psiPS_old = psiPS;
                % Compute new values
                dh = computeDH(hf(2:3),m,mu); % m = 1
                hn(2) = hn(1)*((1/mr)-1i); % m = 1
                Bm = computeBm(m,dh,hn(2)); % m = 1
                psiP = conj(Bm)*P(2);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                [~,~,lF] = evaluateConvergence(psiPS_older,psiP_old,...
                    psiPS_old,psiP,psiPS,thVal,methodFlag);
                while lF
                    m = m + 1;
                    % Store old values
                    psiP_old = psiP;
                    psiPS_older = psiPS_old;
                    psiPS_old = psiPS;
                    % Compute new values
                    hf(1) = hf(2);
                    hf(2) = hf(3);
                    hf(3) = computeHX(hf(1:2),m+1,mu); % m+1 is ok, not m
                    dh = computeDH(hf(2:3),m,mu);
                    hn(3) = computeHX(hn(1:2),m,mr); % m is ok, not m+1
                    Bm = computeBm(m,dh,hn(3));
                    P(3) = computeP(P(1:2),m,x);
                    psiP = conj(Bm)*P(3);
                    psiPS = psiPS + psiP;
                    psiPCS{m+1} = psiPS;
                    P(1) = P(2);
                    P(2) = P(3);
                    hn(1) = hn(2);
                    hn(2) = hn(3);
                    [~,~,lF] = evaluateConvergence(psiPS_older,psiP_old,...
                        psiPS_old,psiP,psiPS,thVal,methodFlag);
                end
                H = -(rho/mu)*exp(-1i*mu*rho)*psiPS;
            end
            % Compute output order, N
            if methodFlag == 4
                psiPCS = cell2mat(psiPCS);
                if thVal ~= inf && thVal ~= -inf
                    cVec = round(diff(psiPCS),thVal);
                else
                    cVec = diff(psiPCS);
                end
                nonNanIndxs = find(~isnan(cVec));
                lastNonZeroChange = find(cVec(nonNanIndxs),1,'last');
                if isempty(lastNonZeroChange)
                    N = 0;
                else
                    N = nonNanIndxs(lastNonZeroChange)+1;
                end
            else
                N = m-1;
            end
            % Compute equivalent thresholds
            for ii = 1:2
                [thVec(ii),~,~] = evaluateConvergence(psiPS_older,...
                    psiP_old,psiPS_old,psiP,psiPS,0,ii); % 0 - dummy value
            end
        otherwise
            error('Invalid METHOD specification.')
    end
end

end

function out = computeHX(hx,m,z)
%COMPUTEHX Spherical-Hankel function of the first kind.
%   Y = COMPUTEHX(A,M,Z) computes the value, Y, of the spherical-hankel
%   function of the first kind of order M and with argument Z. A is a
%   2-element vector containing the values of the function for the same
%   argument but for orders M-2 and M-1, stored in this order.

out = (((2*m-1)/z)*hx(2))-hx(1);

end

function out = computeDH(hf,m,z)
%COMPUTEDH Derivative of the spherical-Hankel function of the first kind.
%   Y = COMPUTEDH(A,M,Z) computes the value, Y, of the derivative of the 
%   spherical-hankel function of the first kind of order M and with 
%   argument Z. A is a 2-element vector containing the values of the 
%   spherical-hankel function of the first kind for the same argument, Z,
%   but for orders M and M+1, stored in this order.

out = -hf(2)+((m/z)*hf(1));

end

function out = computeP(P,m,x)
%COMPUTEP Legendre polynomial
%   Y = COMPUTEP(P,M,X) computes the value, Y, of the Legendre polynomial
%   of degree M and with argument X. P is a 2-element vector containing the 
%   values of the Legendre polynomial for the same argument but for degrees
%   M-2 and M-1, stored in this order.

out = (((2*m)-1)/m)*x*P(2)-(((m-1)/m)*P(1));

end

function out = computeAm(m,hp)

out = ((-1i)^(m-1))*(2*m+1)/hp;

end

function out = computeBm(m,hp,hn)

out = ((2*m)+1)*hn/hp;

end

function [old,new,lF] = evaluateConvergence(psiPS_older,psiP_old,...
    psiPS_old,psiP,psiPS,TH,methodFlag)

switch methodFlag
    case 1
        old = computeConvergenceMetric(psiPS_old,psiPS_older,methodFlag);
        new = computeConvergenceMetric(psiPS,psiPS_old,methodFlag);
        lF = old > TH || new > TH;
    case 2
        old = computeConvergenceMetric(psiP_old,psiPS_old,methodFlag);
        new = computeConvergenceMetric(psiP,psiPS,methodFlag);
        lF = old > TH || new > TH;
    case 3
        old = computeConvergenceMetric(psiPS_old,psiPS_older,methodFlag);
        new = computeConvergenceMetric(psiPS,psiPS_old,methodFlag);
        ro = real(old);
        io = imag(old);
        rn = real(new);
        in = imag(new);
        if pres ~= inf
            ro = round(ro,TH);
            io = round(io,TH);
            rn = round(rn,TH);
            in = round(in,TH);
        end
        lF1 = ro ~= 0 || io ~= 0 || rn ~= 0 || in ~= 0;
        lF2 = ~isnan(ro) && ~isnan(io) && ~isnan(rn) && ~isnan(in);
        lF = lF1 && lF2;
    case 4 % exact
        old = 0; % dummy value
        new = 0; % dummy value
        lF = ~isnan(psiP);
    otherwise
        error('Invalid methodFlag value.')
end

end

function out = computeConvergenceMetric(in1,in2,MF)

switch MF
    case 1 % Cooper and Bauck 1980
        out = abs((in1-in2)/in1);
    case 2 % Duda and Martens 1998
        out = abs(in1)/abs(in2);
    case 3 % Sridhar and Choueiri 2019
        out = in1-in2;
    otherwise
        error('Invalid MF specification.')
end

end
