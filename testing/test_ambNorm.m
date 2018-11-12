% Testing ambisonics normalizations

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
%   Copyright (c) 2018 Princeton University
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

Fs = 96000;
FFTLen = 16384;
freqVec = getFreqVec(Fs,FFTLen);
kLen = 1 + (FFTLen / 2);
kVec = f2k(freqVec(1:kLen));

r = [0.1 0.2 0];
s = [0.2 -0.2 0.1];

Lmax = 16;
Nmax = (Lmax + 1)^2;

ambNorm1 = 'N3D';
ambNorm2 = 'SN3D';
ambNorm3 = 'none';

%% Point source

psi = exp(1i*kVec*norm(r-s))/norm(r-s);
psi1 = zeros(size(kVec));
psi2 = zeros(size(kVec));
psi3 = zeros(size(kVec));
for l = 0:Lmax
    jl = sphericalBesselJ(l,kVec*norm(r));
    hl = sphericalHankelH(l,1,kVec*norm(s));
    coeff1 = 4*pi*1i*kVec.*jl.*hl/ambNormSquared(l,ambNorm1);
    coeff2 = 4*pi*1i*kVec.*jl.*hl/ambNormSquared(l,ambNorm2);
    coeff3 = 4*pi*1i*kVec.*jl.*hl/ambNormSquared(l,ambNorm3);
    for m = -l:l
        Yr1 = ambSphericalHarmonicY(l, m, r, ambNorm1);
        Ys1 = ambSphericalHarmonicY(l, m, s, ambNorm1);
        Yr2 = ambSphericalHarmonicY(l, m, r, ambNorm2);
        Ys2 = ambSphericalHarmonicY(l, m, s, ambNorm2);
        Yr3 = ambSphericalHarmonicY(l, m, r, ambNorm3);
        Ys3 = ambSphericalHarmonicY(l, m, s, ambNorm3);
        psi1 = psi1 + coeff1*Yr1*Ys1;
        psi2 = psi2 + coeff2*Yr2*Ys2;
        psi3 = psi3 + coeff3*Yr3*Ys3;
    end
end

figure()
plot(kVec,abs(psi))
hold all
plot(kVec,abs(psi1))
plot(kVec,abs(psi2),':')
plot(kVec,abs(psi3),'--')
hold off
set(gca,'XScale','log')

%% Plane-wave source

psi = exp(-1i*kVec*dot(s,r));
psi1 = zeros(size(kVec));
psi2 = zeros(size(kVec));
psi3 = zeros(size(kVec));
for l = 0:Lmax
    jl = sphericalBesselJ(l,kVec*norm(r));
    coeff1 = 4*pi*((-1i)^l)*jl/ambNormSquared(l,ambNorm1);
    coeff2 = 4*pi*((-1i)^l)*jl/ambNormSquared(l,ambNorm2);
    coeff3 = 4*pi*((-1i)^l)*jl/ambNormSquared(l,ambNorm3);
    for m = -l:l
        Yr1 = ambSphericalHarmonicY(l, m, r, ambNorm1);
        Ys1 = ambSphericalHarmonicY(l, m, s, ambNorm1);
        Yr2 = ambSphericalHarmonicY(l, m, r, ambNorm2);
        Ys2 = ambSphericalHarmonicY(l, m, s, ambNorm2);
        Yr3 = ambSphericalHarmonicY(l, m, r, ambNorm3);
        Ys3 = ambSphericalHarmonicY(l, m, s, ambNorm3);
        psi1 = psi1 + coeff1*Yr1*Ys1;
        psi2 = psi2 + coeff2*Yr2*Ys2;
        psi3 = psi3 + coeff3*Yr3*Ys3;
    end
end

figure()
plot(kVec,abs(psi))
hold all
plot(kVec,abs(psi1))
plot(kVec,abs(psi2),':')
plot(kVec,abs(psi3),'--')
hold off
set(gca,'XScale','log')

%% Plane-wave decomposition

[vMat,wVec] = loadGridFile('fliege_36');
numPW = length(wVec);
pinvFlag = true; % works for true or false

Lmax = 4;
Nmax = (Lmax + 1)^2;

aIndx = 7;
AMat = zeros(kLen,Nmax);
AMat(:,aIndx) = 1;

muMat1 = a2mu(AMat,vMat,ambNorm1,pinvFlag,wVec);
muMat2 = a2mu(AMat,vMat,ambNorm2,pinvFlag,wVec);
muMat3 = a2mu(AMat,vMat,ambNorm3,pinvFlag,wVec);

AMat1 = mu2a(muMat1,vMat,Lmax,wVec,ambNorm1);
AMat2 = mu2a(muMat2,vMat,Lmax,wVec,ambNorm2);
AMat3 = mu2a(muMat3,vMat,Lmax,wVec,ambNorm3);

figure()
subplot(2,2,1)
plot(kVec,abs(AMat(:,aIndx)))
hold all
plot(kVec,abs(AMat1(:,aIndx)))
plot(kVec,abs(AMat2(:,aIndx)),':')
plot(kVec,abs(AMat3(:,aIndx)),'--')
ylim([0,2])
hold off
set(gca,'XScale','log')
subplot(2,2,2)
plot(kVec,abs(AMat(:,1)))
hold all
plot(kVec,abs(AMat1(:,1)))
plot(kVec,abs(AMat2(:,1)),':')
plot(kVec,abs(AMat3(:,1)),'--')
hold off
ylim([0,2])
set(gca,'XScale','log')
subplot(2,2,3)
plot(kVec,abs(AMat(:,aIndx+2)))
hold all
plot(kVec,abs(AMat1(:,aIndx+2)))
plot(kVec,abs(AMat2(:,aIndx+2)),':')
plot(kVec,abs(AMat3(:,aIndx+2)),'--')
hold off
ylim([0,2])
set(gca,'XScale','log')
subplot(2,2,4)
plot(kVec,abs(AMat(:,Nmax)))
hold all
plot(kVec,abs(AMat1(:,Nmax)))
plot(kVec,abs(AMat2(:,Nmax)),':')
plot(kVec,abs(AMat3(:,Nmax)),'--')
hold off
ylim([0,2])
set(gca,'XScale','log')

%% Plane-wave synthesis - this one only seems to work if Nmax = numPW

[vMat,wVec] = loadGridFile('fliege_144');
numPW = length(wVec);

Lmax = 11;
Nmax = (Lmax + 1)^2; % needs to match numPW

pwIndx = 5;
muMat = zeros(kLen,numPW);
muMat(:,pwIndx) = 1;

AMat1 = mu2a(muMat,vMat,Lmax,wVec,ambNorm1);
AMat2 = mu2a(muMat,vMat,Lmax,wVec,ambNorm2);
AMat3 = mu2a(muMat,vMat,Lmax,wVec,ambNorm3);

muMat1 = a2mu(AMat1,vMat,ambNorm1);
muMat2 = a2mu(AMat2,vMat,ambNorm2);
muMat3 = a2mu(AMat3,vMat,ambNorm3);

figure()
subplot(2,2,1)
plot(kVec,abs(muMat(:,pwIndx)))
hold all
plot(kVec,abs(muMat1(:,pwIndx)))
plot(kVec,abs(muMat2(:,pwIndx)),':')
plot(kVec,abs(muMat3(:,pwIndx)),'--')
ylim([0,2])
hold off
set(gca,'XScale','log')
subplot(2,2,2)
plot(kVec,abs(muMat(:,1)))
hold all
plot(kVec,abs(muMat1(:,1)))
plot(kVec,abs(muMat2(:,1)),':')
plot(kVec,abs(muMat3(:,1)),'--')
hold off
ylim([0,2])
set(gca,'XScale','log')
subplot(2,2,3)
plot(kVec,abs(muMat(:,pwIndx+2)))
hold all
plot(kVec,abs(muMat1(:,pwIndx+2)))
plot(kVec,abs(muMat2(:,pwIndx+2)),':')
plot(kVec,abs(muMat3(:,pwIndx+2)),'--')
hold off
ylim([0,2])
set(gca,'XScale','log')
subplot(2,2,4)
plot(kVec,abs(muMat(:,numPW)))
hold all
plot(kVec,abs(muMat1(:,numPW)))
plot(kVec,abs(muMat2(:,numPW)),':')
plot(kVec,abs(muMat3(:,numPW)),'--')
hold off
ylim([0,2])
set(gca,'XScale','log')

%% Plane-wave synthesis - using pinv for a2mu

[vMat,wVec] = loadGridFile('fliege_36');
numPW = length(wVec);

Lmax = 10; % doesn't work if Lmax is too small
Nmax = (Lmax + 1)^2;

pwIndx = 5;
muMat = zeros(kLen,numPW);
muMat(:,pwIndx) = 1;

AMat1 = mu2a(muMat,vMat,Lmax,wVec,ambNorm1);
AMat2 = mu2a(muMat,vMat,Lmax,wVec,ambNorm2);
AMat3 = mu2a(muMat,vMat,Lmax,wVec,ambNorm3);

muMat1 = a2mu(AMat1,vMat,ambNorm1,true,wVec);
muMat2 = a2mu(AMat2,vMat,ambNorm2,true,wVec);
muMat3 = a2mu(AMat3,vMat,ambNorm3,true,wVec);

figure()
subplot(2,2,1)
plot(kVec,abs(muMat(:,pwIndx)))
hold all
plot(kVec,abs(muMat1(:,pwIndx)))
plot(kVec,abs(muMat2(:,pwIndx)),':')
plot(kVec,abs(muMat3(:,pwIndx)),'--')
ylim([0,2])
hold off
set(gca,'XScale','log')
subplot(2,2,2)
plot(kVec,abs(muMat(:,1)))
hold all
plot(kVec,abs(muMat1(:,1)))
plot(kVec,abs(muMat2(:,1)),':')
plot(kVec,abs(muMat3(:,1)),'--')
hold off
ylim([0,2])
set(gca,'XScale','log')
subplot(2,2,3)
plot(kVec,abs(muMat(:,pwIndx+2)))
hold all
plot(kVec,abs(muMat1(:,pwIndx+2)))
plot(kVec,abs(muMat2(:,pwIndx+2)),':')
plot(kVec,abs(muMat3(:,pwIndx+2)),'--')
hold off
ylim([0,2])
set(gca,'XScale','log')
subplot(2,2,4)
plot(kVec,abs(muMat(:,numPW)))
hold all
plot(kVec,abs(muMat1(:,numPW)))
plot(kVec,abs(muMat2(:,numPW)),':')
plot(kVec,abs(muMat3(:,numPW)),'--')
hold off
ylim([0,2])
set(gca,'XScale','log')

%% Norm-squared calculation

[vMat,wVec] = loadGridFile('fliege_900');

Lmax = 10;
Nmax = (Lmax + 1)^2;

normSq1 = zeros(Nmax,1);
normSq2 = zeros(Nmax,1);
normSq3 = zeros(Nmax,1);
intYSq1 = zeros(Nmax,1);
intYSq2 = zeros(Nmax,1);
intYSq3 = zeros(Nmax,1);
for l = 0:Lmax
    for m = -l:l
        n = getACN(l,m);
        normSq1(n+1,1) = ambNormSquared(l, ambNorm1);
        normSq2(n+1,1) = ambNormSquared(l, ambNorm2);
        normSq3(n+1,1) = ambNormSquared(l, ambNorm3);
        Y1 = ambSphericalHarmonicY(l, m, vMat, ambNorm1);
        Y2 = ambSphericalHarmonicY(l, m, vMat, ambNorm2);
        Y3 = ambSphericalHarmonicY(l, m, vMat, ambNorm3);
        intYSq1(n+1,1) = dot((Y1.*Y1),wVec);
        intYSq2(n+1,1) = dot((Y2.*Y2),wVec);
        intYSq3(n+1,1) = dot((Y3.*Y3),wVec);
    end
end

figure()
subplot(3,1,1)
plot(normSq1)
hold all
plot(intYSq1,'--')
subplot(3,1,2)
plot(normSq2)
hold all
plot(intYSq2,'--')
subplot(3,1,3)
plot(normSq3)
hold all
plot(intYSq3,'--')