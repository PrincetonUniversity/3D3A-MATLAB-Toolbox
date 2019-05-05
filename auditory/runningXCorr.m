function [g,L1,L2,m] = runningXCorr(x1,x2,fS,T)
%RUNNINGXCORR Compute running, normalized cross-correlation.
%   G = RUNNINGXCORR(X1,X2,FS) computes the running cross-correlation of 
%   two input signals, X1 and X2, both sampled at a rate of FS Hz. X1 and 
%   X2 must be vectors of the same length. The output, G, is a matrix of 
%   normalized cross-correlation values at each sample and for each lag, m. 
%   By default, m takes values -M:1:M, where M = floor(FS/1000). The
%   output, G, is normalized to have values between 0 and 1.
%
%   G = RUNNINGXCORR(X1,X2,FS,T) optionally specifies the time constant, in 
%   seconds, of the exponentially-decaying estimation window used to 
%   compute the running cross-correlation. The default value is 0.01. T can
%   take values in the range 1/FS to infinity (specified as inf).
%
%   [G,L1,L2] = RUNNINGXCORR(...) additionally returns autocorrelation 
%   matrices, L1 and L2, corresponding to inputs X1 and X2, respectively.
%
%   [G,L1,L2,M] = RUNNINGXCORR(...) additionally returns the lag vector, M.

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

%   Ref:
%       [1]. Faller and Merimaa (2004) - Source localization in complex 
%       listening situations/ Selection of binaural cues based on 
%       interaural coherence.

narginchk(3,4);

if nargin < 4
    T = 0.01;
end

% Check inputs
validateattributes(x1,{'double'},{'vector','nonempty','nonnan','finite',...
    'real'},'runningXCorr','X1',1);
xLen = length(x1);
validateattributes(x2,{'double'},{'vector','nonempty','nonnan','finite',...
    'real','numel',xLen},'runningXCorr','X2',2);
validateattributes(fS,{'double'},{'scalar','nonempty','nonnan','finite',...
    'real','positive'},'runningXCorr','FS',3);
validateattributes(T,{'double'},{'scalar','nonempty','nonnan','>=',...
    1/fS},'runningXCorr','T',4);

x1 = shiftdim(x1);
x2 = shiftdim(x2);

mB = floor(fS/1000);
m = (-mB:mB).';
mLen = length(m);

% Initialize variables
a12 = zeros(xLen+1,mLen);
a11 = zeros(xLen+1,mLen);
a22 = zeros(xLen+1,mLen);

% Compute running, normalized cross-correlation
alpha = 1/(T*fS);
x1p = [zeros(mB,1); x1; zeros(mB,1)];
x2p = [zeros(mB,1); x2; zeros(mB,1)];
for ii = 1:xLen
    for jj = 1:mLen
        a12(ii+1,jj) = alpha*x1p(ii-max(m(jj),0)+mB)*...
            x2p(ii-max(-m(jj),0)+mB) + (1-alpha)*a12(ii,jj);
        a11(ii+1,jj) = alpha*x1p(ii-max(m(jj),0)+mB)^2 + ...
            (1-alpha)*a11(ii,jj);
        a22(ii+1,jj) = alpha*x2p(ii-max(-m(jj),0)+mB)^2 + ...
            (1-alpha)*a22(ii,jj);
    end
end

L1 = a11(2:end,:);
L2 = a22(2:end,:);
g = a12(2:end,:)./sqrt(L1.*L2 + eps); % Small error term to prevent NaNs

end
