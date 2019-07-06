function varargout = computeSphereHRTFs(S,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

narginchk(1,17);

% Parse and verify inputs
[inputs,defs] = parseInputs(S,varargin);

% Extract parsed inputs
S = inputs.S;
A = inputs.A;
T = inputs.T;
Fs = inputs.Fs;
eL = inputs.eL;
eR = inputs.eR;
R = inputs.R;
causalFlag = inputs.makecausal;
compType = inputs.compType;

% Perform additional checks on inputs
if Fs < 1000
    warning(['Input sampling rate, %d, is low. Note that sampling',...
        ' rate must be specified in Hz.'],Fs)
end
if eL(1) < 0
    error('Invalid az specification for eL. az must be within [0,360).')
end
if eR(1) < 0
    error('Invalid az specification for eR. az must be within [0,360).')
end
if eL(2) > 90
    error('Invalid el specification for eL. el must be within [-90,90].')
end
if eR(2) > 90
    error('Invalid el specification for eR. el must be within [-90,90].')
end
if R <= A
    error('Invalid R specification. R must be > A.')
end

% Initialize variables
irLen = round(T*Fs);
fVec = getFreqVec(Fs,irLen);
nyquistIndx = ceil((irLen+1)/2);
c = getSoundSpeed();
mu = fVec(1:nyquistIndx)*2*pi*A/c;

% Compute angles of incidence from source position data
fprintf('Computing angle(s) of incidence...')
numPos = size(S,1);
[xS,yS,zS] = sofaS2sofaC(S(:,1),S(:,2),ones(numPos,1));
[xEL,yEL,zEL] = sofaS2sofaC(eL(1),eL(2),1);
[xER,yER,zER] = sofaS2sofaC(eR(1),eR(2),1);
thetaL = getCentralAngle([xS,yS,zS],repmat([xEL,yEL,zEL],numPos,1));
thetaR = getCentralAngle([xS,yS,zS],repmat([xER,yER,zER],numPos,1));
fprintf('done!\n')

switch lower(compType{1})
    case 'fixedn'
        switch length(compType)
            case 1
                N = defs{8}{2};
            case 2
                N = compType{2};
            otherwise
                error(['Invalid options specification for ''fixedN'' ',...
                    'computation type.'])
        end
        mVec = (0:N).';
        PL = legendreP(repmat(mVec,1,numPos),...
            repmat(cosd(thetaL.'),N+1,1));
        PR = legendreP(repmat(mVec,1,numPos),...
            repmat(cosd(thetaR.'),N+1,1));
        [hrtfL,hrtfR,NMatL,NMatR,dTL,dTR] = computeFixedNHRTFs(R,A,mu,...
            N,PL,PR,irLen,Fs,c,causalFlag);
    case 'auton1' % Proposed
        switch length(compType)
            case 1
                Pr = defs{8}{2};
            case 2
                Pr = compType{2};
            otherwise
                error(['Invalid options specification for ''fixedN'' ',...
                    'computation type.'])
        end
        
        fprintf(['Starting HRTF computation type auton1 with Pr = ',...
            num2str(Pr),'...\n'])
        [hrtfL,hrtfR,NMatL,NMatR,dTL,dTR] = computeAutoN1HRTFs(R,A,mu,...
            thetaL,thetaR,Pr,irLen,Fs,c,causalFlag);
    case 'auton2' % Duda and Martens
        switch length(compType)
            case 1
                th = defs{8}{2};
            case 2
                th = compType{2};
            otherwise
                error(['Invalid options specification for ''fixedN'' ',...
                    'computation type.'])
        end
        
        [hrtfL,hrtfR,NMatL,NMatR,dTL,dTR] = computeAutoN2HRTFs(R,A,mu,...
            thetaL,thetaR,th,irLen,Fs,c,causalFlag);
    otherwise
        error('Invalid input for computation type.')
end

hL = ifft(hrtfL,irLen,1,'symmetric');
hR = ifft(hrtfR,irLen,1,'symmetric');

% Parse outputs
switch nargout
    case 1
        varargout{1} = hL;
    case 2
        varargout{1} = hL;
        varargout{2} = hR;
    case 3
        varargout{1} = hL;
        varargout{2} = hR;
        varargout{3} = NMatL;
    case 4
        varargout{1} = hL;
        varargout{2} = hR;
        varargout{3} = NMatL;
        varargout{4} = NMatR;
    case 5
        varargout{1} = hL;
        varargout{2} = hR;
        varargout{3} = NMatL;
        varargout{4} = NMatR;
        varargout{5} = dTL;
    case 6
        varargout{1} = hL;
        varargout{2} = hR;
        varargout{3} = NMatL;
        varargout{4} = NMatR;
        varargout{5} = dTL;
        varargout{6} = dTR;
    otherwise
        error('Invalid number of output arguments requested.')
end

end

function [inputs,defs] = parseInputs(S,opts)

p = inputParser;

% Specify defaults
defs = {0.09;0.005;44100;[90,0];[270,0];inf;false;{'fixedN',41}};

% Required inputs
addRequired(p,'S',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite','real','ncols',2,'>=',-90,'<',360},...
    'computeSphereHRTFs','S',1));

% Optional inputs
addParameter(p,'A',defs{1},@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','nonnan','finite','real','positive'},...
    'computeSphereHRTFs','A'));
addParameter(p,'T',defs{2},@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','nonnan','finite','real','positive'},...
    'computeSphereHRTFs','T'));
addParameter(p,'Fs',defs{3},@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','nonnan','finite','real','positive'},...
    'computeSphereHRTFs','Fs'));
addParameter(p,'eL',defs{4},@(x)validateattributes(x,{'double'},...
    {'vector','nonempty','nonnan','finite','real','numel',2,'>=',-90,...
    '<',360},'computeSphereHRTFs','eL'));
addParameter(p,'eR',defs{5},@(x)validateattributes(x,{'double'},...
    {'vector','nonempty','nonnan','finite','real','numel',2,'>=',-90,...
    '<',360},'computeSphereHRTFs','eR'));
addParameter(p,'R',defs{6},@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','nonnan','real','positive'},...
    'computeSphereHRTFs','R'));
addParameter(p,'makecausal',defs{7},@(x)validateattributes(x,...
    {'logical'},{'nonempty','nonnan'},'computeSphereHRTFs',...
    'makecausal, c'));
addParameter(p,'compType',defs{8},@(x)validateattributes(x,{'cell'},...
    {'vector','nonempty'},'computeSphereHRTFs','compType'));

p.CaseSensitive = false;
p.FunctionName = 'computeSphereHRTFs';

parse(p,S,opts{:});

inputs = p.Results;

end

function [hrtfL,hrtfR,NMatL,NMatR,dTL,dTR] = computeFixedNHRTFs(R,A,mu,...
    N,PL,PR,irLen,Fs,c,causalFlag)

nyquistIndx = length(mu);
numPos = size(PL,2);
mVec = (0:N).';
hrtfL = ones(nyquistIndx,numPos);
hrtfR = ones(nyquistIndx,numPos);
NMatL = zeros(nyquistIndx,numPos);
NMatR = zeros(nyquistIndx,numPos);
dTL = zeros(nyquistIndx,numPos);
dTR = zeros(nyquistIndx,numPos);
if R == inf
    for ii = 2:nyquistIndx
        dh = dSphericalHankelH(mVec,1,mu(ii));
        % conj so that negative phase = delay
        psi = diag(conj((2*mVec+1).*((-1i).^(mVec-1))./dh));
        psiPL = psi*PL;
        psiPR = psi*PR;
        NMatL(ii,:) = N*ones(1,numPos);
        NMatR(ii,:) = N*ones(1,numPos);
        psiL = sum(psiPL,1);
        psiR = sum(psiPR,1);
        dTL(ii,:) = psiPL(end,:)./psiL;
        dTR(ii,:) = psiPR(end,:)./psiR;
        if causalFlag
            % extra exp term to make IRs causal
            del = mu(ii)*(c/A)*(irLen/Fs)*0.5;
            hrtfL(ii,:) = (1/(mu(ii)^2))*exp(-1i*del)*psiL;
            hrtfR(ii,:) = (1/(mu(ii)^2))*exp(-1i*del)*psiR;
        else
            hrtfL(ii,:) = (1/(mu(ii)^2))*psiL;
            hrtfR(ii,:) = (1/(mu(ii)^2))*psiR;
        end
    end
else
    rho = R/A;
    for ii = 2:nyquistIndx
        h = sphericalHankelH(mVec,1,mu(ii)*rho);
        dh = dSphericalHankelH(mVec,1,mu(ii));
        psi = diag(conj(((2*mVec)+1).*(h./dh)));
        psiPL = psi*PL;
        psiPR = psi*PR;
        NMatL(ii,:) = N*ones(1,numPos);
        NMatR(ii,:) = N*ones(1,numPos);
        psiL = sum(psiPL,1);
        psiR = sum(psiPR,1);
        dTL(ii,:) = psiPL(end,:)./psiL;
        dTR(ii,:) = psiPR(end,:)./psiR;
        if causalFlag
            hrtfL(ii,:) = (rho/mu(ii))*exp(-1i*mu(ii)*rho)*psiL;
            hrtfR(ii,:) = (rho/mu(ii))*exp(-1i*mu(ii)*rho)*psiR;
        else
            hrtfL(ii,:) = (rho/mu(ii))*psiL(ii,:);
            hrtfR(ii,:) = (rho/mu(ii))*psiR(ii,:);
        end
    end
end

end

function [hrtfL,hrtfR,NMatL,NMatR,dTL,dTR] = computeAutoN1HRTFs(R,A,mu,...
    thetaL,thetaR,Pr,irLen,Fs,c,causalFlag)

nyquistIndx = length(mu);
numPos = size(thetaL,1);
hrtfL = ones(nyquistIndx,numPos);
hrtfR = ones(nyquistIndx,numPos);
NMatL = zeros(nyquistIndx,numPos);
NMatR = zeros(nyquistIndx,numPos);
psiL = zeros(nyquistIndx,numPos);
psiR = zeros(nyquistIndx,numPos);
dTL = zeros(nyquistIndx,numPos);
dTR = zeros(nyquistIndx,numPos);
if R == inf
    clearProgress = repmat('\b',1,6);
    fprintf('Computing HRTFs...%5.1f%%',0);
    for ii = 2:nyquistIndx
        % Progress
        fprintf(clearProgress);
        fprintf('%5.1f%%',(ii-1)*100/nyquistIndx);
        for jj = 1:numPos
            Nmax = 0;
            dh = dSphericalHankelH(Nmax,1,mu(ii));
            % conj so that negative phase = delay
            psi = conj((2*Nmax+1).*((-1i).^(Nmax-1))./dh);
            PL = legendreP(Nmax,cosd(thetaL(jj)));
            psiPL = psi*PL;
            cutoff = 0;
            while ~isnan(psiPL) && cutoff < 2
                oldPsiL = psiL(ii,jj);
                psiL(ii,jj) = psiL(ii,jj) + psiPL;
                tempDTL = psiPL/psiL(ii,jj);
                
                if Pr ~= inf
                    if round(psiL(ii,jj),Pr) == round(oldPsiL,Pr)
                        cutoff = cutoff + 1;
                    else
                        cutoff = 0;
                    end
                else
                    if psiL(ii,jj) == oldPsiL
                        cutoff = cutoff + 1;
                    else
                        cutoff = 0;
                    end
                end
                
                Nmax = Nmax + 1;
                dh = dSphericalHankelH(Nmax,1,mu(ii));
                psi = conj((2*Nmax+1).*((-1i).^(Nmax-1))./dh);
                PL = legendreP(Nmax,cosd(thetaL(jj))); 
                psiPL = psi*PL;
            end
            psiL(ii,jj) = oldPsiL;
            NMatL(ii,jj) = Nmax-2;
            dTL(ii,jj) = tempDTL;
            
            Nmax = 0;
            PR = legendreP(Nmax,cosd(thetaR(jj)));
            psiPR = psi*PR;
            cutoff = 0;
            while ~isnan(psiPR) && cutoff < 2
                oldPsiR = psiR(ii,jj);
                psiR(ii,jj) = psiR(ii,jj) + psiPR;
                tempDTR = psiPR/psiR(ii,jj);
                
                if Pr ~= inf
                    if round(psiR(ii,jj),Pr) == round(oldPsiR,Pr)
                        cutoff = cutoff + 1;
                    else
                        cutoff = 0;
                    end
                else
                    if psiR(ii,jj) == oldPsiR
                        cutoff = cutoff + 1;
                    else
                        cutoff = 0;
                    end
                end
                
                Nmax = Nmax + 1;
                dh = dSphericalHankelH(Nmax,1,mu(ii));
                psi = conj((2*Nmax+1).*((-1i).^(Nmax-1))./dh);
                PR = legendreP(Nmax,cosd(thetaR(jj)));
                psiPR = psi*PR;
            end
            psiR(ii,jj) = oldPsiR;
            NMatR(ii,jj) = Nmax-2;
            dTR(ii,jj) = tempDTR;
        end
        if causalFlag
            % extra exp term to make IRs causal
            del = mu(ii)*(c/A)*(irLen/Fs)*0.5;
            hrtfL(ii,:) = (1/(mu(ii)^2))*exp(-1i*del)*psiL(ii,:);
            hrtfR(ii,:) = (1/(mu(ii)^2))*exp(-1i*del)*psiR(ii,:);
        else
            hrtfL(ii,:) = (1/(mu(ii)^2))*psiL(ii,:);
            hrtfR(ii,:) = (1/(mu(ii)^2))*psiR(ii,:);
        end
    end
    fprintf(clearProgress);
    fprintf('%5.1f%%\n',100);
else
    rho = R/A;
    clearProgress = repmat('\b',1,6);
    fprintf('Computing HRTFs...%5.1f%%',0);
    for ii = 2:nyquistIndx
        % Progress
        fprintf(clearProgress);
        fprintf('%5.1f%%',(ii-1)*100/nyquistIndx);
        for jj = 1:numPos
            Nmax = 0;
            h = sphericalHankelH(Nmax,1,mu(ii)*rho);
            dh = dSphericalHankelH(Nmax,1,mu(ii));
            psi = conj(((2*Nmax)+1).*(h./dh));
            PL = legendreP(Nmax,cosd(thetaL(jj)));
            psiPL = psi*PL;
            cutoff = 0;
            while ~isnan(psiPL) && cutoff < 2
                oldPsiL = psiL(ii,jj);
                psiL(ii,jj) = psiL(ii,jj) + psiPL;
                tempDTL = psiPL/psiL(ii,jj);
                
                if Pr ~= inf
                    if round(psiL(ii,jj),Pr) == round(oldPsiL,Pr)
                        cutoff = cutoff + 1;
                    else
                        cutoff = 0;
                    end
                else
                    if psiL(ii,jj) == oldPsiL
                        cutoff = cutoff + 1;
                    else
                        cutoff = 0;
                    end
                end
                
                Nmax = Nmax + 1;
                h = sphericalHankelH(Nmax,1,mu(ii)*rho);
                dh = dSphericalHankelH(Nmax,1,mu(ii));
                psi = conj(((2*Nmax)+1).*(h./dh));
                PL = legendreP(Nmax,cosd(thetaL(jj)));
                psiPL = psi*PL;
            end
            psiL(ii,jj) = oldPsiL;
            NMatL(ii,jj) = Nmax-2;
            dTL(ii,jj) = tempDTL;
            
            Nmax = 0;
            PR = legendreP(Nmax,cosd(thetaR(jj)));
            psiPR = psi*PR;
            cutoff = 0;
            while ~isnan(psiPR) && cutoff < 2
                oldPsiR = psiR(ii,jj);
                psiR(ii,jj) = psiR(ii,jj) + psiPR;
                tempDTR = psiPR/psiR(ii,jj);
                
                if Pr ~= inf
                    if round(psiR(ii,jj),Pr) == round(oldPsiR,Pr)
                        cutoff = cutoff + 1;
                    else
                        cutoff = 0;
                    end
                else
                    if psiR(ii,jj) == oldPsiR
                        cutoff = cutoff + 1;
                    else
                        cutoff = 0;
                    end
                end
                
                Nmax = Nmax + 1;
                h = sphericalHankelH(Nmax,1,mu(ii)*rho);
                dh = dSphericalHankelH(Nmax,1,mu(ii));
                psi = conj(((2*Nmax)+1).*(h./dh));
                PR = legendreP(Nmax,cosd(thetaR(jj)));
                psiPR = psi*PR;
            end
            psiR(ii,jj) = oldPsiR;
            NMatR(ii,jj) = Nmax-2;
            dTR(ii,jj) = tempDTR;
        end
        if causalFlag
            hrtfL(ii,:) = (rho/mu(ii))*exp(-1i*mu(ii)*rho)*psiL(ii,:);
            hrtfR(ii,:) = (rho/mu(ii))*exp(-1i*mu(ii)*rho)*psiR(ii,:);
        else
            hrtfL(ii,:) = (rho/mu(ii))*psiL(ii,:);
            hrtfR(ii,:) = (rho/mu(ii))*psiR(ii,:);
        end
    end
    fprintf(clearProgress);
    fprintf('%5.1f%%\n',100);
end

end

function [hrtfL,hrtfR,NMatL,NMatR,dTL,dTR] = computeAutoN2HRTFs(R,A,mu,...
    thetaL,thetaR,th,irLen,Fs,c,causalFlag)

nyquistIndx = length(mu);
numPos = size(thetaL,1);
hrtfL = ones(nyquistIndx,numPos);
hrtfR = ones(nyquistIndx,numPos);
NMatL = zeros(nyquistIndx,numPos);
NMatR = zeros(nyquistIndx,numPos);
psiL = zeros(nyquistIndx,numPos);
psiR = zeros(nyquistIndx,numPos);
dTL = zeros(nyquistIndx,numPos);
dTR = zeros(nyquistIndx,numPos);
if R == inf    
    for ii = 2:nyquistIndx
        for jj = 1:numPos
            Nmax = 0;
            dh = dSphericalHankelH(Nmax,1,mu(ii));
            % conj so that negative phase = delay
            psi = conj((2*Nmax+1).*((-1i).^(Nmax-1))./dh);
            PL = legendreP(Nmax,cosd(thetaL(jj)));
            psiPL = psi*PL;
            psiL(ii,jj) = psiL(ii,jj) + psiPL;
            oldDTL = 1;
            dTL(ii,jj) = psiPL/psiL(ii,jj);
            while oldDTL > th && dTL(ii,jj) > th                
                Nmax = Nmax + 1;
                dh = dSphericalHankelH(Nmax,1,mu(ii));
                psi = conj((2*Nmax+1).*((-1i).^(Nmax-1))./dh);
                PL = legendreP(Nmax,cosd(thetaL(jj)));
                psiPL = psi*PL;
                psiL(ii,jj) = psiL(ii,jj) + psiPL;
                oldDTL = dTL(ii,jj);
                dTL(ii,jj) = psiPL/psiL(ii,jj);
            end
            NMatL(ii,jj) = Nmax;
            
            Nmax = 0;
            PR = legendreP(Nmax,cosd(thetaR(jj)));
            psiPR = psi*PR;
            psiR(ii,jj) = psiR(ii,jj) + psiPR;
            oldDTR = 1;
            dTR(ii,jj) = psiPR/psiR(ii,jj);
            while oldDTR > th && dTR(ii,jj) > th
                Nmax = Nmax + 1;
                dh = dSphericalHankelH(Nmax,1,mu(ii));
                psi = conj((2*Nmax+1).*((-1i).^(Nmax-1))./dh);
                PR = legendreP(Nmax,cosd(thetaR(jj)));
                psiPR = psi*PR;
                psiR(ii,jj) = psiR(ii,jj) + psiPR;
                oldDTR = dTR(ii,jj);
                dTR(ii,jj) = psiPR/psiR(ii,jj);
            end
            NMatR(ii,jj) = Nmax;
        end
        if causalFlag
            % extra exp term to make IRs causal
            del = mu(ii)*(c/A)*(irLen/Fs)*0.5;
            hrtfL(ii,:) = (1/(mu(ii)^2))*exp(-1i*del)*psiL(ii,:);
            hrtfR(ii,:) = (1/(mu(ii)^2))*exp(-1i*del)*psiR(ii,:);
        else
            hrtfL(ii,:) = (1/(mu(ii)^2))*psiL(ii,:);
            hrtfR(ii,:) = (1/(mu(ii)^2))*psiR(ii,:);
        end
    end
else
    rho = R/A;
    for ii = 2:nyquistIndx
        for jj = 1:numPos
            Nmax = 0;
            h = sphericalHankelH(Nmax,1,mu(ii)*rho);
            dh = dSphericalHankelH(Nmax,1,mu(ii));
            psi = conj(((2*Nmax)+1).*(h./dh));
            PL = legendreP(Nmax,cosd(thetaL(jj)));
            psiPL = psi*PL;
            psiL(ii,jj) = psiL(ii,jj) + psiPL;
            oldDTL = 1;
            dTL(ii,jj) = psiPL/psiL(ii,jj);
            while oldDTL > th && dTL(ii,jj) > th
                Nmax = Nmax + 1;
                h = sphericalHankelH(Nmax,1,mu(ii)*rho);
                dh = dSphericalHankelH(Nmax,1,mu(ii));
                psi = conj(((2*Nmax)+1).*(h./dh));
                PL = legendreP(Nmax,cosd(thetaL(jj)));
                psiPL = psi*PL;
                psiL(ii,jj) = psiL(ii,jj) + psiPL;
                oldDTL = dTL(ii,jj);
                dTL(ii,jj) = psiPL/psiL(ii,jj);
            end
            NMatL(ii,jj) = Nmax;
            
            Nmax = 0;
            PR = legendreP(Nmax,cosd(thetaR(jj)));
            psiPR = psi*PR;
            psiR(ii,jj) = psiR(ii,jj) + psiPR;
            oldDTR = 1;
            dTR(ii,jj) = psiPR/psiR(ii,jj);
            while oldDTR > th && dTR(ii,jj) > th
                Nmax = Nmax + 1;
                h = sphericalHankelH(Nmax,1,mu(ii)*rho);
                dh = dSphericalHankelH(Nmax,1,mu(ii));
                psi = conj(((2*Nmax)+1).*(h./dh));
                PR = legendreP(Nmax,cosd(thetaR(jj)));
                psiPR = psi*PR;
                psiR(ii,jj) = psiR(ii,jj) + psiPR;
                oldDTR = dTR(ii,jj);
                dTR(ii,jj) = psiPR/psiR(ii,jj);
            end
            NMatR(ii,jj) = Nmax;
        end
        if causalFlag
            hrtfL(ii,:) = (rho/mu(ii))*exp(-1i*mu(ii)*rho)*psiL(ii,:);
            hrtfR(ii,:) = (rho/mu(ii))*exp(-1i*mu(ii)*rho)*psiR(ii,:);
        else
            hrtfL(ii,:) = (rho/mu(ii))*psiL(ii,:);
            hrtfR(ii,:) = (rho/mu(ii))*psiR(ii,:);
        end
    end
end

end
