function A = ambPointSource(L,s0,k,r,ambNorm,hpfFlag)
%AMBPOINTSOURCE Ambisonics potentials for an omnidirectional point source.
%   A = AMBPOINTSOURCE(L,S0,K,R) computes the ambisonics potentials A, for
%   angular wavenumber K and up to order L, given a point-source signal
%   originating from S0 and evaluated at R (both given in Cartesian
%   coordinates).
%
%   A = AMBPOINTSOURCE(L,S0,K,R,AMBNORM) additionally specifies the
%   ambisonics normalization convention to use. By default, N3D is used.
%
%   A = AMBPOINTSOURCE(L,S0,K,R,AMBNORM,HPFFLAG) additionally specifies the
%   near-field compensation high-pass filter to use. By default, no filter
%   is applied. The options are: true, false, 'fixed', or 'em32'.
%
%   See also AMBRADIALFILTER, AMBPLANEWAVE.

%   ==============================================================================
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
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   ==============================================================================

narginchk(4,6);

% Uses N3D by default
if nargin < 5 || isempty(ambNorm)
    ambNorm = 'N3D';
end

% No HPF by default
if nargin < 6 || isempty(hpfFlag)
    hpfFlag = false;
end

d = s0 - r; % Effective source position
kLen = length(k);
distFilters = zeros(kLen,L+1);
for ll = 0:L
    distFilters(:,ll+1) = ambRadialFilter(ll,k,norm(d));
end

% Apply order-dependent high-pass filters
if hpfFlag
    if ischar(hpfFlag) && strcmpi(hpfFlag,'em32')
        if L > 4
            warning('EM32 near-field HPFs are undefined for Li > 4.');
        end
        kRef = interp1(2:4,f2k([400 1000 1800]),0:L,'nearest','extrap'); % basically a lookup table
        for ll = 2:L
            magResp = 1 - 1./sqrt(1+(k/kRef(ll+1)).^ll);
            distFilters(:,ll+1) = distFilters(:,ll+1).*magResp;
        end
    else
        if ischar(hpfFlag) && strcmpi(hpfFlag,'fixed')
            baseFreq = 200; % (200 * l) Hz
        else
            baseFreq = 40/norm(d); % (40 * l / d) Hz
        end
        
        for ll = 1:L
            kRef = f2k(baseFreq)*ll;
            magResp = 1 - 1./sqrt(1+(k/kRef).^ll);
            distFilters(:,ll+1) = distFilters(:,ll+1).*magResp;
        end
    end
end

Y = getAmbYMatrix(d, L, ambNorm);
N = (L + 1)^2;
A = zeros(kLen,N);
for nn = 0:(N-1)
    [ll,~] = getAmbOrder(nn);
    A(:,nn+1) = Y(nn+1) * distFilters(:,ll+1);
end

end