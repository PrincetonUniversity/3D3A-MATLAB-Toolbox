function Qa = variableYawAmbRotationMatrix(alpha, maxOrder)
%VARIABLEYAWAMBROTATIONMATRIX Ambisonics rotation in yaw.
%   Q = VARIABLEYAWAMBROTATIONMATRIX(A, L) computes the ambisonic rotation
%   coefficients matrix Q, up to ambisonics order L, for a rotation of A 
%   radians yaw.

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

Qa = zeros(HOATerms);
for ii = 1:HOATerms
    for jj = 1:HOATerms
        if (nList(ii) == nList(jj)) && (abs(mList(ii)) == abs(mList(jj)))
            if mList(ii)*mList(jj) >= 0
                Qa(ii,jj) = cos(mList(jj)*alpha);
            else
                Qa(ii,jj) = sin(mList(jj)*alpha);
            end
        end
    end
end

end