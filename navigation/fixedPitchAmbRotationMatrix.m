function Qb = fixedPitchAmbRotationMatrix(maxOrder)
%FIXEDPITCHAMBROTATIONMATRIX Ambisonics rotation of 90 degrees pitch.
%   Q = FIXEDPITCHAMBROTATIONMATRIX(L) computes the ambisonic rotation
%   coefficients matrix Q, up to ambisonics order L, for a rotation of 90
%   degrees pitch.

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
%   Copyright (c) 2017 Princeton University
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

%   References:
%     [1] Gumerov and Duraiswami (2005) Fast Multipole Methods for the
%         Helmholtz Equation in Three Dimensions.
%     [2] Zotter (2009) Analysis and Synthesis of Sound-Radiation with
%         Spherical Arrays.

HOATerms = (maxOrder + 1)^2;
[nList, mList] = getAmbOrder(0:HOATerms-1);

Q = zeros((2*maxOrder + 1)^2);

% Step 1
for n = 0:(2*maxOrder)
    Pn = legendre(n,0);
    for s = -n:n
        % Eq. 189 [2], m = 0
        Q(getACN(n,0)+1,getACN(n,s)+1) = ((-1)^s) * sqrt(2-(~s)) ...
            * sqrt(factorial(n-abs(s))/factorial(n+abs(s))) * Pn(abs(s)+1);
    end
end

% Step 2
for m = 1:maxOrder
    for n = m:(2*maxOrder - m)
        coeff = sqrt(2-(~m))/(2*sphCoefficientB(n+1,m-1)*sqrt(2-(~(m-1))));
        for s = m:n
            % Eq. 190 [2]
            temp = (sphCoefficientB(n+1,s-1)/sqrt(2-(~(s-1))))*Q(getACN(n+1,m-1)+1,getACN(n+1,s-1)+1) ...
                - (sphCoefficientB(n+1,-s-1)/sqrt(2-(~(s+1))))*Q(getACN(n+1,m-1)+1,getACN(n+1,s+1)+1);
            Q(getACN(n, m)+1,getACN(n,s)+1) = coeff*(sqrt(2-(~s))*temp ...
                + 2*sphCoefficientA(n,s)*Q(getACN(n+1,m-1)+1,getACN(n+1,s)+1));
        end
    end
end

% Step 3
for n = 1:maxOrder
    for m = 1:n
        for s = 0:(m-1)
            % Eq. 191 [2]
            Q(getACN(n,m)+1,getACN(n,s)+1) = ((-1)^(m+s))*Q(getACN(n,s)+1,getACN(n,m)+1);
        end
    end
end

% Step 4
for n = 1:maxOrder
    for m = 1:n
        for s = 1:n
            Q(getACN(n,-m)+1,getACN(n,-s)+1) = Q(getACN(n,m)+1,getACN(n,s)+1);
        end
    end
end

% Step 5
kTable = zeros(HOATerms);
for ii = 1:HOATerms
    for jj = 1:HOATerms
        if nList(ii) == nList(jj)
            if (mList(ii) >= 0) && (mList(jj) >= 0)
                kTable(ii,jj) = ~mod(nList(ii) + mList(ii) + mList(jj), 2);
            elseif (mList(ii) < 0) && (mList(jj) < 0)
                kTable(ii,jj) = ~mod(nList(ii) + mList(ii) + mList(jj) + 1, 2);
            end
        end
    end
end

Qb = kTable.*Q(1:HOATerms,1:HOATerms);

end