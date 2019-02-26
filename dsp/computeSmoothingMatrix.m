function M = computeSmoothingMatrix(FFTLen, FRAC, METHOD, WINTYPE)
%COMPUTESMOOTHINGMATRIX Fractional-octave smoothing matrix.
%   M = COMPUTESMOOTHINGMATRIX(HLEN,N,METHOD,WINTYPE) computes a smoothing
%   matrix M for smoothing a transfer function of length HLEN with
%   1/N-octave smoothing, using the specified smoothing METHOD and WINTYPE.
%
%   The returned matrix M will be FLOOR(1+HLEN/2)-by-HLEN.
%
%   See also FRACTIONALOCTAVESMOOTH.

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

%   References:
%     [1] Hatziantoniou and Mourjopoulos (2000) Generalized
%         Fractional-Octave Smoothing of Audio and Acoustic Responses.
%     [2] Tylka et al. (2017) A Generalized Method for Fractional-Octave
%         Smoothing of Transfer Functions that Preserves Log-Frequency
%         Symmetry.

specLen = floor(1 + FFTLen/2);

if FRAC == 0
    warning('No smoothing applied.')
    M = eye(specLen,FFTLen);
    return;
end

switch(lower(METHOD))
    case 'tylka'
        M = smoothingMatrix_tylka(FFTLen, FRAC, WINTYPE);
    case {'hatz','hatziantoniou'}
        M = smoothingMatrix_hatz(FFTLen, FRAC, WINTYPE);
    otherwise
        M = eye(specLen,FFTLen);
end

end

function M = smoothingMatrix_tylka(FFTLen, FRAC, WINTYPE)

specLen = floor(1 + FFTLen/2);
M = zeros(specLen, FFTLen);

M(1,1) = 1;
for n = 1:(specLen-2)
    fL = n*2^(-1/2/FRAC);
    nL = floor(fL);
    fH = n*2^(1/2/FRAC);
    nH = ceil(fH);
    
    % Upper limit correction
    if nH > (specLen - 1)
        nH = specLen - 1;
        nL = floor((n^2)/nH);
    end
    if fH > (FFTLen/2)
        fH = FFTLen/2;
        fL = (n^2)/fH;
    end
    
    % Window function computation
    m = nL:nH;
    winLen = nH - nL + 1;
    tempWinVec = zeros(winLen,1);
    if winLen == 1
        winVec = 1;
    else
        mL = clip(m-0.5, [fL, fH]); % Lower limit of integral
        mH = clip(m+0.5, [fL, fH]); % Upper limit of integral
        mask = (mL - mH) ~= 0;
        switch(lower(WINTYPE))
            case {'rectangular','rect'}
                tempWinVec = (log2(mH/n) - log2(mL/n));
            case {'hanning','hann'}
                term1 = log2(mH/n) - log2(mL/n);
                term2 = sin(2*pi*FRAC*log2(mH/n)) - sin(2*pi*FRAC*log2(mL/n));
                tempWinVec = term1/2 + term2/(4*pi*FRAC);
        end
        tempWinVec = tempWinVec.*mask;
        winVec = tempWinVec/sum(tempWinVec);
    end
    M(n+1,m+1) = winVec;
end

M(specLen, specLen) = 1;

end

function M = smoothingMatrix_hatz(FFTLen, FRAC, WINTYPE)

specLen = floor(1 + FFTLen/2);
M = zeros(specLen, FFTLen);

Q = 1 / (2^(0.5 / FRAC) - 2^(-0.5 / FRAC));

M(1,1) = 1;
for n = 1:(specLen-1)
    % Window half-width calculation
    m = floor((0.5 * n) / Q);
    if m > n
        m = n;
    end
    
    % Window function computation
    indx = (n-m):(n+m);
    winLen = 2*m + 1;
    if winLen == 1
        winVec = 1;
    else
        switch(lower(WINTYPE))
            case {'rectangular','rect'}
                tempWinVec = ones(1,winLen);
            case {'hanning','hann'}
                tempWinVec = hann(winLen);
        end
        winVec = tempWinVec/sum(tempWinVec);
    end
    M(n+1,indx+1) = winVec;
end

end