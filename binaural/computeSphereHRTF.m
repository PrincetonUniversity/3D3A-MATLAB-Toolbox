function [H,N,thVec] = computeSphereHRTF(a,r,theta,f,METHOD)
%COMPUTESPHEREHRTF Spherical-head HRTFs.
%   Detailed explanation goes here

narginchk(4,5);

if nargin < 5
    METHOD{1} = 'cooperbauck1980';
    METHOD{2} = 1e-3; % threshold
end

c = getSoundSpeed();
mu = (2*pi*f*a)/c;
thVec = zeros(2,1);
if mu == 0
    H = 1;
    N = 0;
else
    x = cosd(theta);
    rho = r/a;
    switch lower(METHOD{1})
        case 'fixedn'
            thVal = 0; % Needed to compute equivalent thresholds in the end
            N = METHOD{2}; % order

            % Compute first 2 terms
            hf = zeros(3,1); % Initialize hf to store 3 most recent terms
            P = zeros(3,1); % Initialize P to store 3 most recent terms
            hf(1) = exp(1i*mu)/(1i*mu); % m = 0
            hf(2) = hf(1)*((1/mu)-1i); % m = 1
            hf(3) = computeH(hf(1:2),2,mu); % m = 2
            P(1) = 1; % m = 0
            P(2) = x; % m = 1
            
            % Initialize more variables
            psiPCS = zeros(N+1,1);
            psiPS_old = 0;
            psiPS = 0;
            if rho == inf
                m = 0;
                dh = computeDH(hf(1:2),m,mu); % m = 0
                Am = computeAm(m,dh); % m = 0
                psiP = conj(Am)*P(1);
                psiPS = psiPS + psiP;
                psiPCS(m+1) = psiPS;
                m = m + 1;
                psiP_old = psiP;
                psiPS_older = psiPS_old;
                psiPS_old = psiPS;
                dh = computeDH(hf(2:3),m,mu); % m = 1
                Am = computeAm(m,dh); % m = 1
                psiP = conj(Am)*P(2);
                psiPS = psiPS + psiP;
                psiPCS(m+1) = psiPS;
                for m = 2:N 
                    hf(1) = hf(2);
                    hf(2) = hf(3);
                    hf(3) = computeH(hf(1:2),m+1,mu);
                    psiP_old = psiP;
                    psiPS_older = psiPS_old;
                    psiPS_old = psiPS;
                    dh = computeDH(hf(2:3),m,mu);
                    Am = computeAm(m,dh);
                    P(3) = computeP(P(1:2),m,x);
                    psiP = conj(Am)*P(3);
                    psiPS = psiPS + psiP;
                    psiPCS(m+1) = psiPS;
                    P(1) = P(2);
                    P(2) = P(3);
                end
                H = (1/mu^2)*psiPS;
            else
                mr = mu*rho;
                hn = zeros(3,1);
                m = 0;
                dh = computeDH(hf(1:2),m,mu); % m = 0
                hn(1) = exp(1i*mr)/(1i*mr); % m = 0
                Bm = computeBm(m,dh,hn(1)); % m = 0
                psiP = conj(Bm)*P(1);
                psiPS = psiPS + psiP;
                psiPCS(m+1) = psiPS;
                m = m + 1;
                psiP_old = psiP;
                psiPS_older = psiPS_old;
                psiPS_old = psiPS;
                dh = computeDH(hf(2:3),m,mu); % m = 1
                hn(2) = hn(1)*((1/mr)-1i); % m = 1
                Bm = computeBm(m,dh,hn(2)); % m = 1
                psiP = conj(Bm)*P(2);
                psiPS = psiPS + psiP;
                psiPCS(m+1) = psiPS;
                for m = 2:N
                    hf(1) = hf(2);
                    hf(2) = hf(3);
                    hf(3) = computeH(hf(1:2),m+1,mu);
                    psiP_old = psiP;
                    psiPS_older = psiPS_old;
                    psiPS_old = psiPS;
                    dh = computeDH(hf(2:3),m,mu);
                    hn(3) = computeH(hn(1:2),m,mr);
                    Bm = computeBm(m,dh,hn(3));
                    P(3) = computeP(P(1:2),m,x);
                    psiP = conj(Bm)*P(3);
                    psiPS = psiPS + psiP;
                    psiPCS(m+1) = psiPS;
                    P(1) = P(2);
                    P(2) = P(3);
                    hn(1) = hn(2);
                    hn(2) = hn(3);
                end
                H = -(rho/mu)*exp(-1i*mu*rho)*psiPS;
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
            thVal = METHOD{2}; % threshold
            % Compute first 2 terms
            hf = zeros(3,1); % Initialize hf to store 3 most recent terms
            P = zeros(3,1); % Initialize P to store 3 most recent terms
            hf(1) = exp(1i*mu)/(1i*mu); % m = 0
            hf(2) = hf(1)*((1/mu)-1i); % m = 1
            hf(3) = computeH(hf(1:2),2,mu); % m = 2
            P(1) = 1; % m = 0
            P(2) = x; % m = 1
            
            % Initialize more variables
            psiPS_old = 0;
            psiPS = 0;
            if rho == inf
                m = 0;
                dh = computeDH(hf(1:2),m,mu); % m = 0
                Am = computeAm(m,dh); % m = 0
                psiP = conj(Am)*P(1);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                m = m + 1;
                psiP_old = psiP;
                psiPS_older = psiPS_old;
                psiPS_old = psiPS;
                dh = computeDH(hf(2:3),m,mu); % m = 1
                Am = computeAm(m,dh); % m = 1
                psiP = conj(Am)*P(2);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                [~,~,lF] = computeMultiThresh(psiPS_older,psiP_old,...
                    psiPS_old,psiP,psiPS,thVal,methodFlag);
                while lF
                    m = m + 1;
                    hf(1) = hf(2);
                    hf(2) = hf(3);
                    hf(3) = computeH(hf(1:2),m+1,mu);
                    psiP_old = psiP;
                    psiPS_older = psiPS_old;
                    psiPS_old = psiPS;
                    dh = computeDH(hf(2:3),m,mu);
                    Am = computeAm(m,dh);
                    P(3) = computeP(P(1:2),m,x);
                    psiP = conj(Am)*P(3);
                    psiPS = psiPS + psiP;
                    psiPCS{m+1} = psiPS;
                    P(1) = P(2);
                    P(2) = P(3);
                    [~,~,lF] = computeMultiThresh(psiPS_older,...
                        psiP_old,psiPS_old,psiP,psiPS,thVal,methodFlag);
                end
                H = (1/mu^2)*psiPS;
            else
                mr = mu*rho;
                hn = zeros(3,1);
                m = 0;
                dh = computeDH(hf(1:2),m,mu); % m = 0
                hn(1) = exp(1i*mr)/(1i*mr); % m = 0
                Bm = computeBm(m,dh,hn(1)); % m = 0
                psiP = conj(Bm)*P(1);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                m = m + 1;
                psiP_old = psiP;
                psiPS_older = psiPS_old;
                psiPS_old = psiPS;
                dh = computeDH(hf(2:3),m,mu); % m = 1
                hn(2) = hn(1)*((1/mr)-1i); % m = 1
                Bm = computeBm(m,dh,hn(2)); % m = 1
                psiP = conj(Bm)*P(2);
                psiPS = psiPS + psiP;
                psiPCS{m+1} = psiPS;
                [~,~,lF] = computeMultiThresh(psiPS_older,psiP_old,...
                    psiPS_old,psiP,psiPS,thVal,methodFlag);
                while lF
                    m = m + 1;
                    hf(1) = hf(2);
                    hf(2) = hf(3);
                    hf(3) = computeH(hf(1:2),m+1,mu);
                    psiP_old = psiP;
                    psiPS_older = psiPS_old;
                    psiPS_old = psiPS;
                    dh = computeDH(hf(2:3),m,mu);
                    hn(3) = computeH(hn(1:2),m,mr);
                    Bm = computeBm(m,dh,hn(3));
                    P(3) = computeP(P(1:2),m,x);
                    psiP = conj(Bm)*P(3);
                    psiPS = psiPS + psiP;
                    psiPCS{m+1} = psiPS;
                    P(1) = P(2);
                    P(2) = P(3);
                    hn(1) = hn(2);
                    hn(2) = hn(3);
                    [~,~,lF] = computeMultiThresh(psiPS_older,psiP_old,...
                        psiPS_old,psiP,psiPS,thVal,methodFlag);
                end
                H = -(rho/mu)*exp(-1i*mu*rho)*psiPS;
            end
            if methodFlag == 4
                psiPCS = cell2mat(psiPCS);
                pres = log10(1/thVal);
                if pres ~= inf && pres ~= -inf
                    cVec = round(diff(psiPCS),pres);
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
        otherwise
            error('Invalid METHOD specification.')
    end
    % Compute equivalent thresholds
    for ii = 1:2
        [ro,rn,~] = computeMultiThresh(psiPS_older,psiP_old,...
            psiPS_old,psiP,psiPS,thVal,ii);
        thVec(ii) = max([ro,rn]);
    end
end

end

function out = computeDH(hf,m,z)

out = -hf(2)+((m/z)*hf(1));

end

function out = computeAm(m,hp)

out = ((-1i)^(m-1))*(2*m+1)/hp;

end

function out = computeH(hf,m,z)

out = (((2*m-1)/z)*hf(2))-hf(1);

end

function [old,new,lF] = computeMultiThresh(psiPS_older,psiP_old,...
    psiPS_old,psiP,psiPS,TH,methodFlag)

switch methodFlag
    case 1
        old = computeThresh(psiPS_old,psiPS_older,methodFlag);
        new = computeThresh(psiPS,psiPS_old,methodFlag);
        lF = old > TH || new > TH;
    case 2
        old = computeThresh(psiP_old,psiPS_old,methodFlag);
        new = computeThresh(psiP,psiPS,methodFlag);
        lF = old > TH || new > TH;
    case 3
        old = computeThresh(psiPS_old,psiPS_older,methodFlag);
        new = computeThresh(psiPS,psiPS_old,methodFlag);
        pres = log10(1/TH);
        ro = real(old);
        io = imag(old);
        rn = real(new);
        in = imag(new);
        if pres ~= inf
            ro = round(ro,pres);
            io = round(io,pres);
            rn = round(rn,pres);
            in = round(in,pres);
        end
        lF1 = ro ~= 0 || io ~= 0 || rn ~= 0 || in ~= 0;
        lF2 = ~isnan(ro) && ~isnan(io) && ~isnan(rn) && ~isnan(in);
        lF = lF1 && lF2;
    case 4 % exact
        old = 0; % dummy
        new = 0; % dummy
        lF = ~isnan(psiP);
    otherwise
        error('Invalid methodFlag value.')
end

end

function out = computeP(P,m,x)

out = (((2*m)-1)/m)*x*P(2)-(((m-1)/m)*P(1));

end

function out = computeBm(m,hp,hn)

out = ((2*m)+1)*hn/hp;

end

function out = computeThresh(in1,in2,MF)

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
