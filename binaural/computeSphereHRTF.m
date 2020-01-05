function [H,N,thVec] = computeSphereHRTF(a,r,theta,f,varargin)
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
%   the METHOD to use to compute the HRTF. METHOD must be a cell array and
%   can take the following values:
%       1. {'sridharchoueiri2019',P,NI} where P can be any finite integer  
%       or inf. The value of precision, P, refers to the number of decimal 
%       places to use when determining the order, N, of the calculation, 
%       with P = inf corresponding to maximum precision. NI specifies the 
%       number of partial sums that are checked against a threshold when 
%       determining the order, N. NI must be a positive integer. If METHOD
%       is not specified, {'sridharchoueiri2019',inf,2} is used.
%
%       2. {'fixedn',N} where N is the order and must be a non-negative
%       integer. In this case, the second output, N, will be the same as 
%       this specified input order.
%
%       3. {'cooperbauck1980',TH,NI} where TH can be any finite, real
%       number, and refers to the threshold for determining the order, N, 
%       of the calculation. In the Fortran code published by Cooper and 
%       Bauck [2], TH is set to 0.001. In this case, the specified value of 
%       TH will also be returned as the first element of the output vector, 
%       T. NI specifies the number of partial sums that are checked against 
%       TH when determining the order, N. NI must be a positive integer.
%
%       4. {'dudamartens1998',TH,NI} where TH is as defined in 3, above. In
%       this case, the specified value of TH will also be returned as the
%       second element of the output vector, T. NI specifies the number of 
%       partial sums that are checked against TH when determining the 
%       order, N. NI must be a positive integer.
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
%   [H,N,T] = COMPUTESPHEREHRTF(A,R,THETA,F,METHOD,NORMLOC) optionally
%   allows the normalization used for computing the HRTF to be specified.
%   The two options are:
%       1. 'center' (default) - HRTF is computed by normalizing the
%       scattered pressure by the free-field pressure computed at a
%       position corresponding to the location of the center of the sphere.
%
%       2. 'ear' - HRTF is computed by normalizing the scattered pressure
%       by the free-field pressure computed at the appropriate ear
%       position.
%   To specify NORMLOC with default METHOD, specify METHOD as {}.
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

narginchk(4,6);

if nargin < 6
    NORMLOC = 'center';
else
    NORMLOC = varargin{2};
end

defaultMETHOD = {'sridharchoueiri2019',inf,2};
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
    'finite','real','nonnegative','<',360},'computeSphereHRTF','THETA',3);
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
end
validateattributes(NORMLOC,{'char'},{'scalartext','nonempty'},...
    'computeSphereHRTF','NORMLOC',6);

% Initialize common variables
c = getSoundSpeed();
mu = (2*pi*f*a)/c;
thVec = zeros(2,1);
x = cosd(theta);
rho = r/a;

% Main calculation
if mu == 0
    if rho == inf
        H = 1;
        N = 0;
        thVec = [inf;inf]; % Threshold values don't matter in this case
    else
        P = zeros(3,1); % Initialize P to store 3 most recent terms
        
        % Initialize sums
        psiPS = 0;
        switch lower(METHOD{1})
            case 'fixedn'
                numChkTerms = 1;
                psiPSVec = zeros(numChkTerms+1,1);
                psiPVec = zeros(numChkTerms+1,1);
                
                N = METHOD{2}; % Extract order from input
                m = 0;
                P(3) = 1; % Specify most recent P value (i.e., for m = 0)
                Bm = ((2*m)+1)/((m+1)*rho^m); % m = 0
                psiP = Bm*P(3);
                psiPS = psiPS + psiP;
                termIndx = 2;
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
                    Bm = ((2*m)+1)/((m+1)*rho^m);
                    psiP = Bm*P(3);
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
                
                % Store final HRTF value
                psiPVec = circshift(psiPVec,1);
                psiPSVec = circshift(psiPSVec,1);
                H = psiPSVec(numChkTerms+1);
            case {'cooperbauck1980','dudamartens1998',...
                    'sridharchoueiri2019','exact','maxn'}
                switch lower(METHOD{1})
                    case 'cooperbauck1980'
                        methodFlag = 1;
                    case 'dudamartens1998'
                        methodFlag = 2;
                    case {'sridharchoueiri2019','exact'} % Only for mu = 0
                        methodFlag = 3;
                    case 'maxn' % Only for mu = 0
                        methodFlag = 3;
                        METHOD{2} = inf;
                    otherwise
                        error('Invalid METHOD specification.')
                end
                psiPSVec = zeros(numChkTerms+1,1);
                psiPVec = zeros(numChkTerms+1,1);
                
                thVal = METHOD{2}; % Extract threshold from input
                m = 0;
                P(3) = 1; % Specify most recent P value (i.e., for m = 0)
                Bm = ((2*m)+1)/((m+1)*rho^m); % m = 0
                psiP = Bm*P(3);
                psiPS = psiPS + psiP;
                termIndx = 2;
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
                while lF
                    m = m + 1;
                    Bm = ((2*m)+1)/((m+1)*rho^m);
                    psiP = Bm*P(3);
                    psiPS = psiPS + psiP;
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
                
                % Compute output order, N
                N = m-termIndx+1;
                
                % Store final HRTF value
                psiPVec = circshift(psiPVec,1);
                psiPSVec = circshift(psiPSVec,1);
                H = psiPSVec(1);
            otherwise
                error('Invalid METHOD specification.')
        end
        
        % Compute equivalent thresholds
        for ii = 1:2
            thVec(ii) = computeEquivThreshold(psiPSVec,psiPVec,ii);
        end
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
                termIndx = 2;
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
                psiPVec = circshift(psiPVec,1);
                psiPSVec = circshift(psiPSVec,1);
                switch lower(NORMLOC)
                    case 'center'
                        H = (1/mu^2)*psiPSVec(numChkTerms+1);
                    case 'ear'
                        H = (1/mu^2)*exp(-1i*mu*x)*psiPSVec(numChkTerms+1);
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
                termIndx = 2;
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
                psiPVec = circshift(psiPVec,1);
                psiPSVec = circshift(psiPSVec,1);
                switch lower(NORMLOC)
                    case 'center'
                        H = -(rho/mu)*psiPSVec(numChkTerms+1);
                    case 'ear'
                        H = -(rho/mu)*exp(1i*mu*(sqrt(1+rho^2-2*rho*...
                            x)))*psiPSVec(numChkTerms+1);
                    otherwise
                        error('Invalid NORMLOC specification.')
                end
            end
        case {'cooperbauck1980','dudamartens1998','sridharchoueiri2019',...
                'exact','maxn'}
            switch lower(METHOD{1})
                case 'cooperbauck1980'
                    methodFlag = 1;
                case 'dudamartens1998'
                    methodFlag = 2;
                case 'sridharchoueiri2019'
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
                termIndx = 2;
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
                while lF
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
                psiPVec = circshift(psiPVec,1);
                psiPSVec = circshift(psiPSVec,1);
                switch lower(NORMLOC)
                    case 'center'
                        H = (1/mu^2)*psiPSVec(1);
                    case 'ear'
                        H = (1/mu^2)*exp(-1i*mu*x)*psiPSVec(1);
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
                termIndx = 2;
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
                while lF
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
                psiPVec = circshift(psiPVec,1);
                psiPSVec = circshift(psiPSVec,1);
                switch lower(NORMLOC)
                    case 'center'
                        H = -(rho/mu)*psiPSVec(1);
                    case 'ear'
                        H = -(rho/mu)*exp(1i*mu*(sqrt(1+rho^2-2*rho*...
                            x)))*psiPSVec(1);
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
                    N = nonNanIndxs(lastNonZeroChange)+1;
                end
                [H,~,thVec] = computeSphereHRTF(a,r,theta,f,...
                    {'fixedn',N},NORMLOC);
            else
                N = m-termIndx+1;
                
                % Compute equivalent thresholds
                for ii = 1:2
                    thVec(ii) = computeEquivThreshold(psiPSVec,psiPVec,ii);
                end
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
%   LF = EVALUATECONVERGENCE(PSV,PV,TH,MF) returns a loog flag used to
%   indicate if the series solution has converged. This is don by using the
%   input parameters PSV and/or PV, and a convergence threshold/numerical
%   precision value, TH.

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
%COMPUTECONVERGENCEMETRIC Convergence metric for determing order.
%   Y = COMPUTECONVERGENCEMETRIC(X1,X2,FLAG) computes the value, Y, of the
%   convergence metric corresponding to the method identified by the flag,
%   MF, given input parameters X1 and X2.

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

function out = computeEquivThreshold(psiPSVec,psiPVec,methodFlag)

switch methodFlag
    case 1
        out = computeConvergenceMetric(psiPSVec(2),psiPSVec(1),methodFlag);
    case 2
        out = computeConvergenceMetric(psiPVec(2),psiPSVec(2),methodFlag);
    otherwise
        error('Invalid methodFlag value.')
end

end
