function [H,N,thVec] = computeSphereHRTF(a,r,theta,f,varargin)
%COMPUTESPHEREHRTF Analytically-computed HRTF for a sphere.
%   [H,N,T] = COMPUTESPHEREHRTF(A,R,THETA,F) analytically computes the HRTF
%   value, H, for a sphere of radius A (in meters), source at distance R
%   (in meters) with angle of incidence THETA (in degrees), and for
%   frequency F (in Hz) using the method described by Sridhar and Choueiri
%   [1] with a numerical precision value of inf (i.e., maximum precision).
%   All inputs must be real-valued scalars. THETA can take values in the 
%   range [0,180]. Also returned are the computed order, N, used to compute 
%   H, and a 2-by-1 vector, T, of corresponding threshold estimates for use 
%   with the HRTF computation methods described by Cooper and Bauck [2], 
%   and Duda and Martens [3], in that order. The GETCENTRALANGLE function 
%   may be useful for computing THETA.
%
%   [H,N,T] = COMPUTESPHEREHRTF(A,R,THETA,F,METHOD) optionally specifies
%   the METHOD to use to compute the HRTF. METHOD must be a cell array and
%   can take the following values:
%       1. {'sridharchoueiri2020',P,NI} where P can be any finite integer  
%       or inf. The value of precision, P, refers to the number of decimal 
%       places to use when determining the order, N, of the calculation, 
%       with P = inf corresponding to maximum precision. NI specifies the 
%       number of partial sums that are checked against a threshold when 
%       determining the order, N. NI must be a positive integer. If METHOD
%       is not specified, {'sridharchoueiri2020',inf,4} is used.
%
%       2. {'fixedn',N} where N is the order and must be a non-negative
%       integer. In this case, the second output, N, will be the same as 
%       this specified input order.
%
%       3. {'cooperbauck1980',TH,NI} where TH can be any finite, real
%       number, and refers to the threshold for determining the order, N, 
%       of the calculation. In the Fortran code published by Cooper and 
%       Bauck [2], TH is set to 0.001. NI specifies the number of partial 
%       sums that are checked against TH when determining the order, N. 
%       NI must be a positive integer.
%
%       4. {'dudamartens1998',TH,NI} where TH and NI are as defined in 3, 
%       above.
%
%       5. {'exact',P} where P is as defined in 1, above, computes the
%       smallest order, N, that guarantees a numerically exact (to within
%       precision, P) value of H.
%
%       6. {'maxn'} computes the maximum permissible value of order, N, for
%       the given set of input parameters. The returned value of H is exact
%       but this setting is different from {'exact',inf} in that whereas
%       the latter identifies the smallest value of N that gives a 
%       numerically exact solution, {'maxn'} identifies the largest value 
%       of N that does so, noting that a larger N will produce NaNs. Note
%       that {'maxn'} gives the same result as {'exact',inf} only when F =
%       0.
%
%       7. {'formulan'} computes N using the approximation formulas
%       provided by Sridhar and Choueiri [1] and then uses the 'fixedn'
%       approach defined in 2, above.
%
%   [H,N,T] = COMPUTESPHEREHRTF(A,R,THETA,F,METHOD,NORMLOC) optionally
%   allows the normalization used for computing the HRTF to be specified.
%   The three options are:
%       1. 'center' (default) - HRTF is computed by normalizing the
%       scattered pressure by the free-field pressure computed at a
%       position corresponding to the location of the center of the sphere.
%
%       2. 'ear' - HRTF is computed by normalizing the scattered pressure
%       by the free-field pressure computed at the appropriate ear
%       position.
%
%       3. 'center_tdc' - same as 'center' except for time delay 
%       compensation (tdc) based on the appropriate ear position.
%       
%       4. 'center_tdc_corr',FMAX - same as 'center_tdc' except the 
%       time-delay compensation accuracy is improved by using creeping wave 
%       theory. For more, see the work by Sridhar and Choueiri [1]. For 
%       this option, the frequency, FMAX, at which the creeping wave 
%       assumption is taken to be valid should be specified in Hz. 
%       Typically, when computing spherical-head HRTFs as discrete-time 
%       filters, FMAX may be specified as the Nyquist frequency for the 
%       chosen sampling rate.
%   To specify NORMLOC with default METHOD, specify METHOD as {}.
%
%   [H,N,T] = COMPUTESPHEREHRTF(A,R,THETA,F,METHOD,NORMLOC,'built-in')
%   optionally allows MATLAB's built-in functions to be used to compute the
%   values of special functions. The built-in functions employ Miller's
%   recurrence algorithms to avoid rounding errors in the calculation of
%   some of the special functions. Note: Specifying this option can result
%   in very long computation times.
%
%   See also GETCENTRALANGLE, COMPUTESPHEREITD, COMPUTESPHEREILD.

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

% Refs:
%   [1]. Sridhar and Choueiri (2020) - A Formula for Quickly Computing the 
%   RS-HRTF.
%   [2]. Cooper and Bauck (1980) - On Acoustical Specification of Natural 
%   Stereo Imaging.
%   [3]. Duda and Martens (1998) - Range dependence of the response of a 
%   spherical head model.

narginchk(4,8);

global funcFlag
if any(strcmpi(varargin,'built-in'))
    funcFlag = true;
else
    funcFlag = false;
end

if nargin < 6
    NORMLOC = 'center';
else
    NORMLOC = varargin{2};
end

defaultMETHOD = {'sridharchoueiri2020',inf,4};
if nargin < 5
    METHOD = defaultMETHOD;
else
    METHOD = varargin{1};
end

% Check inputs
validateattributes(a,{'numeric'},{'scalar','nonempty','nonnan','finite',...
    'positive','real'},'computeSphereHRTF','A',1);
validateattributes(r,{'numeric'},{'scalar','nonempty','nonnan'},...
    'computeSphereHRTF','R',2);
validateattributes(theta,{'numeric'},{'scalar','nonempty','nonnan',...
    'finite','real','nonnegative','<=',180},'computeSphereHRTF','THETA',3);
validateattributes(f,{'numeric'},{'scalar','nonempty','nonnan','finite',...
    'nonnegative','real'},'computeSphereHRTF','F',4);
validateattributes(METHOD,{'cell'},{'2d'},'computeSphereHRTF','METHOD',5);
if isempty(METHOD)
    METHOD = defaultMETHOD;
else
    validateattributes(METHOD{1},{'char'},{'scalartext','nonempty'},...
        'computeSphereHRTF','first element of METHOD',5);
    methodLen = length(METHOD);
    if methodLen > 1
        validateattributes(METHOD{2},{'numeric'},{'scalar','nonempty',...
            'nonnan'},'computeSphereHRTF','second element of METHOD',5);
    end
    
    if methodLen > 2
        validateattributes(METHOD{3},{'numeric'},{'scalar','nonempty',...
            'nonnan','integer','positive'},'computeSphereHRTF',...
            'third element of METHOD',5);
        numChkTerms = METHOD{3};
    else
        numChkTerms = 2;
    end
    
    % The following is provided for backwards compatibility.
    if strcmpi(METHOD{1},'sridharchoueiri2019')
        METHOD{1} = 'sridharchoueiri2020';
    end
end
validateattributes(NORMLOC,{'char'},{'scalartext','nonempty'},...
    'computeSphereHRTF','NORMLOC',6);

% Initialize common variables
c = getSoundSpeed();
mu = (2*pi*f*a)/c;
thVec = zeros(2,1);
x = cosd(theta);
rho = r/a;
g = sqrt(rho^2-(2*rho*x)+1);
N_abort = 500; % Protection to prevent infinite loops
if strcmpi(NORMLOC,'center_tdc_corr')
    if length(varargin) < 3
        error(['FMAX must be specified for NORMLOC option',...
            ' ''center_tdc_corr''']);
    else
        FMAX = varargin{3};
        validateattributes(FMAX,{'numeric'},{'scalar','real','nonempty',...
            'nonnan','positive'},'computeSphereHRTF','FMAX',7);
        mu_cw = 2*pi*(FMAX)*a/c;
        corFac = 1/(1+0.5094*(2*mu_cw^2)^(-1/3));
    end
end

% Main calculation
if mu == 0
    if rho == inf
        H = 1;
        N = 0;
        thVec = [inf;inf]; % Threshold values don't matter in this case
    else
        if theta == 0
            H = (2*rho/g)-(rho*log(rho/(rho-1)));
        else
            H = (2*rho/g)-(rho*log((1-(rho*x)+g)/(rho*(1-x))));
        end
        N = 0;
        thVec = [inf;inf]; % Threshold values don't matter in this case
    end
else % mu > 0
    % Perform preliminary calculations
    % Note: hf is the variable storing the value of the spherical-Hankel 
    % function (of the first kind) with mu as the argument. P is the 
    % variable storing the value of the Legendre polynomial.
    
    hf = zeros(3,1); % Initialize hf to store 3 most recent terms
    P = zeros(3,1); % Initialize P to store 3 most recent terms
    hf(2) = exp(1i*mu)/(1i*mu); % m = 0
    hf(3) = hf(2)*((1/mu)-1i); % m = 1
    
    % Initialize sums
    psiPS = 0;
    switch lower(METHOD{1})
        case 'fixedn'
            numChkTerms = 1;
            % numChkTerms + 1 is needed to compute equivalent thresholds
            psiPSVec = zeros(numChkTerms+1,1);
            psiPVec = zeros(numChkTerms+1,1);
            
            N = METHOD{2}; % extract order from input
            if rho == inf
                m = 0;
                P(3) = 1; % Specify most recent P value (i.e., for m = 0)
                dh = computeDH(hf(2:3),m,mu); % m = 0
                Am = computeAm(m,dh); % m = 0
                psiP = conj(Am)*P(3);
                psiPS = psiPS + psiP;
                termIndx = 1;
                psiPVec(termIndx) = psiP;
                psiPSVec(termIndx) = psiPS;
                if termIndx == (numChkTerms+1)
                    termIndx = termIndx - 1;
                    psiPVec = circshift(psiPVec,-1);
                    psiPSVec = circshift(psiPSVec,-1);  
                end
                P = circshift(P,-1); % Move current value out of P(3)
                P(3) = x; % Update current value (i.e., for m = 1)
                for m = 1:N
                    hf = circshift(hf,-1);
                    hf(3) = computeHX(hf(1:2),m+1,mu); % m+1 is ok, not m
                    dh = computeDH(hf(2:3),m,mu);
                    Am = computeAm(m,dh);
                    psiP = conj(Am)*P(3);
                    psiPS = psiPS + psiP;
                    termIndx = termIndx + 1;
                    psiPVec(termIndx) = psiP;
                    psiPSVec(termIndx) = psiPS;
                    if termIndx == (numChkTerms+1)
                        termIndx = termIndx - 1;
                        psiPVec = circshift(psiPVec,-1);
                        psiPSVec = circshift(psiPSVec,-1);
                    end
                    P = circshift(P,-1); % Move current value out of P(3)
                    P(3) = computeP(P(1:2),m+1,x); % Update current value
                end
                psiPVec = circshift(psiPVec,1); % Undo last circshift
                psiPSVec = circshift(psiPSVec,1); % Undo last circshift
                switch lower(NORMLOC)
                    case 'center'
                        H = (1/mu^2)*psiPSVec(numChkTerms+1);
                    case 'ear'
                        H = (1/mu^2)*exp(-1i*mu*x)*psiPSVec(numChkTerms+1);
                    % ear_to for backwards compatibility.
                    case {'ear_to','center_tdc','center_tdc_corr'} 
                        if theta <= 90
                            delay = mu*x;
                        else
                            if strcmpi(NORMLOC,'center_tdc_corr')
                                delay = -mu*deg2rad(theta-90)/corFac;
                            else
                                delay = -mu*deg2rad(theta-90);
                            end
                        end
                        H = (1/mu^2)*exp(-1i*delay)*...
                            psiPSVec(numChkTerms+1);
                    otherwise
                        error('Invalid NORMLOC specification.')
                end
            else
                mr = mu*rho;
                hn = zeros(3,1); % Init. hn to store 3 most recent terms
                % Note: hn is the variable storing the value of the 
                % spherical-Hankel function (of the first kind) with mr as 
                % the argument.
                hn(3) = exp(1i*mr)/(1i*mr); % m = 0
                
                m = 0;
                P(3) = 1; % Specify most recent P value (i.e., for m = 0)
                dh = computeDH(hf(2:3),m,mu); % m = 0
                Bm = computeBm(m,dh,hn(3)); % m = 0
                psiP = conj(Bm)*P(3);
                psiPS = psiPS + psiP;
                termIndx = 1;
                psiPVec(termIndx) = psiP;
                psiPSVec(termIndx) = psiPS;
                if termIndx == (numChkTerms+1)
                    termIndx = termIndx - 1;
                    psiPVec = circshift(psiPVec,-1);
                    psiPSVec = circshift(psiPSVec,-1);  
                end
                P = circshift(P,-1); % Move current value out of P(3)
                P(3) = x; % Update current value (i.e., for m = 1)
                hn = circshift(hn,-1);
                hn(3) = hn(2)*((1/mr)-1i); % m = 1
                for m = 1:N
                    hf = circshift(hf,-1);
                    hf(3) = computeHX(hf(1:2),m+1,mu); % m+1 is ok, not m
                    dh = computeDH(hf(2:3),m,mu);
                    Bm = computeBm(m,dh,hn(3));
                    psiP = conj(Bm)*P(3);
                    psiPS = psiPS + psiP;
                    termIndx = termIndx + 1;
                    psiPVec(termIndx) = psiP;
                    psiPSVec(termIndx) = psiPS;
                    if termIndx == (numChkTerms+1)
                        termIndx = termIndx - 1;
                        psiPVec = circshift(psiPVec,-1);
                        psiPSVec = circshift(psiPSVec,-1);
                    end
                    P = circshift(P,-1); % Move current value out of P(3)
                    P(3) = computeP(P(1:2),m+1,x); % Update current value
                    hn = circshift(hn,-1);
                    hn(3) = computeHX(hn(1:2),m+1,mr);
                end
                psiPVec = circshift(psiPVec,1); % Undo last circshift
                psiPSVec = circshift(psiPSVec,1); % Undo last circshift
                switch lower(NORMLOC)
                    case 'center'
                        H = -(rho/mu)*exp(1i*mu*rho)*...
                            psiPSVec(numChkTerms+1);
                    case 'ear'                      
                        H = -(g/(rho*mu))*exp(1i*mu*g)*...
                            psiPSVec(numChkTerms+1);
                    % ear_to for backwards compatibility.
                    case {'ear_to','center_tdc','center_tdc_corr'}
                        theta0 = acosd(1/rho);
                        hatG = sqrt(rho^2-1);
                        if theta < theta0
                            delay = mu*g;
                        else
                            if strcmpi(NORMLOC,'center_tdc_corr')
                                delay = mu*hatG + ...
                                    mu*deg2rad(theta-theta0)/corFac;
                            else
                                delay = mu*(hatG + deg2rad(theta-theta0));
                            end 
                        end
                        H = -(rho/mu)*exp(1i*delay)*...
                            psiPSVec(numChkTerms+1);
                    otherwise
                        error('Invalid NORMLOC specification.')
                end
            end
        case {'cooperbauck1980','dudamartens1998','sridharchoueiri2020',...
                'exact','maxn'}
            switch lower(METHOD{1})
                case 'cooperbauck1980'
                    methodFlag = 1;
                case 'dudamartens1998'
                    methodFlag = 2;
                case 'sridharchoueiri2020'
                    methodFlag = 3;
                case 'exact'
                    methodFlag = 4;
                    numChkTerms = 1;
                case 'maxn'
                    methodFlag = 5;
                    METHOD{2} = inf;
                    numChkTerms = 1;
                otherwise
                    error('Invalid METHOD specification.')
            end
            psiPSVec = zeros(numChkTerms+1,1);
            psiPVec = zeros(numChkTerms+1,1);
            
            thVal = METHOD{2}; % Extract threshold from input
            if rho == inf
                m = 0;
                P(3) = 1; % Specify most recent P value (i.e., for m = 0)
                dh = computeDH(hf(2:3),m,mu); % m = 0
                Am = computeAm(m,dh); % m = 0
                psiP = conj(Am)*P(3);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                termIndx = 1;
                psiPVec(termIndx) = psiP;
                psiPSVec(termIndx) = psiPS;
                if termIndx == (numChkTerms+1)
                    lF = evaluateConvergence(psiPSVec,psiPVec,thVal,...
                        methodFlag);
                    termIndx = termIndx - 1;
                    psiPVec = circshift(psiPVec,-1);
                    psiPSVec = circshift(psiPSVec,-1);  
                else
                    lF = true;
                end
                P = circshift(P,-1); % Move current value out of P(3)
                P(3) = x; % Update current value (i.e., for m = 1)
                while lF && (m < N_abort)
                    m = m + 1;
                    hf = circshift(hf,-1);
                    hf(3) = computeHX(hf(1:2),m+1,mu); % m+1 is ok, not m
                    dh = computeDH(hf(2:3),m,mu);
                    Am = computeAm(m,dh);
                    psiP = conj(Am)*P(3);
                    psiPS = psiPS + psiP;
                    psiPCS{m+1} = psiPS;
                    termIndx = termIndx + 1;
                    psiPVec(termIndx) = psiP;
                    psiPSVec(termIndx) = psiPS;
                    if termIndx == (numChkTerms+1)
                        lF = evaluateConvergence(psiPSVec,psiPVec,thVal,...
                            methodFlag);
                        termIndx = termIndx - 1;
                        psiPVec = circshift(psiPVec,-1);
                        psiPSVec = circshift(psiPSVec,-1);
                    else
                        lF = true;
                    end
                    P = circshift(P,-1); % Move current value out of P(3)
                    P(3) = computeP(P(1:2),m+1,x); % Update current value
                end
                psiPVec = circshift(psiPVec,1); % Undo last circshift
                psiPSVec = circshift(psiPSVec,1); % Undo last circshift
                switch lower(NORMLOC)
                    case 'center'
                        H = (1/mu^2)*psiPSVec(1);
                    case 'ear'
                        H = (1/mu^2)*exp(-1i*mu*x)*psiPSVec(1);
                    % ear_to for backwards compatibility.
                    case {'ear_to','center_tdc','center_tdc_corr'}
                        if theta <= 90
                            delay = mu*x;
                        else
                            if strcmpi(NORMLOC,'center_tdc_corr')
                                delay = -mu*deg2rad(theta-90)/corFac;
                            else
                                delay = -mu*deg2rad(theta-90);
                            end
                        end
                        H = (1/mu^2)*exp(-1i*delay)*psiPSVec(1);
                    otherwise
                        error('Invalid NORMLOC specification.')
                end
            else
                mr = mu*rho;
                hn = zeros(3,1); % Init. hn to store 3 most recent terms
                % Note: hn is the variable storing the value of the 
                % spherical-Hankel function (of the first kind) with mr as 
                % the argument.
                hn(3) = exp(1i*mr)/(1i*mr); % m = 0
                
                m = 0;
                P(3) = 1; % Specify most recent P value (i.e., for m = 0)
                dh = computeDH(hf(2:3),m,mu); % m = 0
                Bm = computeBm(m,dh,hn(3)); % m = 0
                psiP = conj(Bm)*P(3);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                termIndx = 1;
                psiPVec(termIndx) = psiP;
                psiPSVec(termIndx) = psiPS;
                if termIndx == (numChkTerms+1)
                    lF = evaluateConvergence(psiPSVec,psiPVec,thVal,...
                        methodFlag);
                    termIndx = termIndx - 1;
                    psiPVec = circshift(psiPVec,-1);
                    psiPSVec = circshift(psiPSVec,-1);  
                else
                    lF = true;
                end
                P = circshift(P,-1); % Move current value out of P(3)
                P(3) = x; % Update current value (i.e., for m = 1)
                hn = circshift(hn,-1);
                hn(3) = hn(2)*((1/mr)-1i); % m = 1
                while lF && (m < N_abort)
                    m = m + 1;
                    hf = circshift(hf,-1);
                    hf(3) = computeHX(hf(1:2),m+1,mu); % m+1 is ok, not m
                    dh = computeDH(hf(2:3),m,mu);
                    Bm = computeBm(m,dh,hn(3));
                    psiP = conj(Bm)*P(3);
                    psiPS = psiPS + psiP;
                    psiPCS{m+1} = psiPS;
                    termIndx = termIndx + 1;
                    psiPVec(termIndx) = psiP;
                    psiPSVec(termIndx) = psiPS;
                    if termIndx == (numChkTerms+1)
                        lF = evaluateConvergence(psiPSVec,psiPVec,thVal,...
                            methodFlag);
                        termIndx = termIndx - 1;
                        psiPVec = circshift(psiPVec,-1);
                        psiPSVec = circshift(psiPSVec,-1);
                    else
                        lF = true;
                    end
                    P = circshift(P,-1); % Move current value out of P(3)
                    P(3) = computeP(P(1:2),m+1,x); % Update current value
                    hn = circshift(hn,-1);
                    hn(3) = computeHX(hn(1:2),m+1,mr);
                end
                psiPVec = circshift(psiPVec,1); % Undo last circshift
                psiPSVec = circshift(psiPSVec,1); % Undo last circshift
                switch lower(NORMLOC)
                    case 'center'
                        H = -(rho/mu)*exp(1i*mu*rho)*psiPSVec(1);
                    case 'ear'                      
                        H = -(g/(rho*mu))*exp(1i*mu*g)*psiPSVec(1);
                    % ear_to for backwards compatibility.
                    case {'ear_to','center_tdc','center_tdc_corr'}
                        theta0 = acosd(1/rho);
                        hatG = sqrt(rho^2-1);
                        if theta < theta0
                            delay = mu*g;
                        else
                            if strcmpi(NORMLOC,'center_tdc_corr')
                                delay = mu*hatG + ...
                                    mu*deg2rad(theta-theta0)/corFac;
                            else
                                delay = mu*(hatG + deg2rad(theta-theta0));
                            end
                        end
                        H = -(rho/mu)*exp(1i*delay)*psiPSVec(1);
                    otherwise
                        error('Invalid NORMLOC specification.')
                end
            end
            
            % Compute output order, N and equivalent thresholds
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
                    N = nonNanIndxs(lastNonZeroChange);
                end
                [H,~,thVec] = computeSphereHRTF(a,r,theta,f,{'fixedn',...
                    N},NORMLOC);
            else
                N = m-termIndx;
            end
            
            if m == N_abort
                warning(['Hard limit of %d iterations reached. ',...
                    'Iteration stopped.'],N_abort)
            end
        case 'formulan'
            % The formula is applicable for rho >= 2 only. For rho < 2, the
            % formula with rho = 2 may be used to obtain satisfactory 
            % results.
            if rho < 2
                rho = 2;
            end
            
            if rho == inf
                alphaVal = 1.41;
                betaVal = 2.73;
                gammaVal = 0.8;
            else
                alphaVal = (1.41*rho+3.9)/(rho-1.36);
                betaVal = (2.73*rho-4.75)/(rho-1.21);
                gammaVal = (0.8*rho-1.01)/(rho-1.45);
            end
            N = round(alphaVal+(betaVal*mu^gammaVal));
            [H,~,thVec] = computeSphereHRTF(a,r,theta,f,{'fixedn',N},...
                NORMLOC);
        otherwise
            error('Invalid METHOD specification.')
    end
    
    % Compute equivalent thresholds
    if ~strcmpi(METHOD{1},'formulan')
        for ii = 1:2
            thVec(ii) = computeEquivThreshold(psiPSVec,psiPVec,ii);
        end
    end
end

end

function out = computeHX(hx,m,z)
%COMPUTEHX Spherical-Hankel function of the first kind.
%   Y = COMPUTEHX(A,M,Z) computes the value, Y, of the spherical-hankel
%   function of the first kind of order M and with argument Z. A is a
%   2-element vector containing the values of the function for the same
%   argument but for orders M-2 and M-1, stored in this order.

global funcFlag
if funcFlag
    out = sphericalHankelH(m,1,z);
else
    out = (((2*m-1)/z)*hx(2))-hx(1);
end

end

function out = computeDH(hf,m,z)
%COMPUTEDH Derivative of the spherical-Hankel function of the first kind.
%   Y = COMPUTEDH(A,M,Z) computes the value, Y, of the derivative of the 
%   spherical-hankel function of the first kind of order M and with 
%   argument Z. A is a 2-element vector containing the values of the 
%   spherical-hankel function of the first kind for the same argument, Z,
%   but for orders M and M+1, stored in this order.

global funcFlag
if funcFlag
    out = dSphericalHankelH(m,1,z);
else
    out = -hf(2)+((m/z)*hf(1));
end

end

function out = computeP(P,m,x)
%COMPUTEP Legendre polynomial
%   Y = COMPUTEP(P,M,X) computes the value, Y, of the Legendre polynomial
%   of degree M and with argument X. P is a 2-element vector containing the 
%   values of the Legendre polynomial for the same argument but for degrees
%   M-2 and M-1, stored in this order.

global funcFlag
if funcFlag
    out = legendreP(m,x);
else
    out = (((2*m)-1)/m)*x*P(2)-(((m-1)/m)*P(1));
end

end

function out = computeAm(m,hp)
%COMPUTEAM Factor in summation term in solution corresponding to infinitely 
%distant source. See the paper by Sridhar and Choueiri [1].

out = ((-1i)^(m-1))*(2*m+1)/hp;

end

function out = computeBm(m,hp,hn)
%COMPUTEBM Factor in summation term in solution corresponding to source at 
%a finite distance. See the paper by Sridhar and Choueiri [1].

out = ((2*m)+1)*hn/hp;

end

function lF = evaluateConvergence(psiPSVec,psiPVec,TH,methodFlag)
%EVALUATECONVERGENCE Determine if series solution has converged.
%   LF = EVALUATECONVERGENCE(PSV,PV,TH,MF) returns a loop flag, LF, used to
%   indicate if the series solution has converged. This is done by using 
%   the input parameters PSV and/or PV, and a convergence threshold/
%   numerical precision value, TH.

psiPVecLen = length(psiPVec);
psiPSVecLen = length(psiPSVec);

switch methodFlag
    case 1
        metricVec = zeros(psiPSVecLen-1,1);
        for ii = 1:(psiPSVecLen-1)
            metricVec(ii) = computeConvergenceMetric(psiPSVec(ii+1),...
                psiPSVec(ii),methodFlag);
        end
        lF = any(metricVec > TH);
    case 2
        if psiPVecLen ~= psiPSVecLen
            error('PSV and PV must have the same length when MF = 2.');
        else
            metricVec = zeros(psiPSVecLen-1,1);
            for ii = 1:(psiPSVecLen-1)
                metricVec(ii) = computeConvergenceMetric(psiPVec(ii+1),...
                    psiPSVec(ii+1),methodFlag);
            end
            lF = any(metricVec > TH);
        end
    case 3
        metricVec = zeros(psiPSVecLen-1,1);
        for ii = 1:(psiPSVecLen-1)
            metricVec(ii) = computeConvergenceMetric(psiPSVec(ii+1),...
                psiPSVec(ii),methodFlag);
        end
        realVec = real(metricVec);
        imagVec = imag(metricVec);
        if TH ~= inf
            realVec = round(realVec,TH);
            imagVec = round(imagVec,TH);
        end
        lF1 = any([realVec;imagVec] ~= 0);
        lF2 = ~any(isnan([realVec;imagVec]));
        lF = lF1 && lF2;
    case {4,5} % exact/maxn
        lF = ~isnan(psiPVec(psiPVecLen));
    otherwise
        error('Invalid methodFlag value.')
end

end

function out = computeConvergenceMetric(in1,in2,MF)
%COMPUTECONVERGENCEMETRIC Convergence metric for determing truncation 
%order.
%   Y = COMPUTECONVERGENCEMETRIC(X1,X2,MF) computes the value, Y, of the
%   convergence metric corresponding to the method identified by the flag,
%   MF, given input parameters X1 and X2.

switch MF
    case 1 % Cooper and Bauck 1980
        out = abs((in1-in2)/in1);
    case 2 % Duda and Martens 1998
        out = abs(in1)/abs(in2);
    case 3 % Sridhar and Choueiri 2020
        out = in1-in2;
    otherwise
        error('Invalid MF specification.')
end

end

function out = computeEquivThreshold(psiPSVec,psiPVec,methodFlag)
%COMPUTEEQUIVTHRESHOLD Estimate equivalent threshold value.
%   O = COMPUTEEQUIVTHRESHOLD(PSV,PV,MF) returns an estimate of the 
%   equivalent threshold value that would need to be used by either the 
%   cooperbauck1980 (MF = 1) or the dudamartens1998 (MF = 2) methods in 
%   order to obtain the PSV and PV vectors provided as inputs.

switch methodFlag
    case 1
        out = computeConvergenceMetric(psiPSVec(2),psiPSVec(1),methodFlag);
    case 2
        out = computeConvergenceMetric(psiPVec(2),psiPSVec(2),methodFlag);
    otherwise
        error('Invalid methodFlag value.')
end

end
