function [rV, rE] = gerzon1992(R, G)
%GERZON1992 Gerzon's localization vectors.
%   [RV, RE] = GERZON1992(R, G) computes Gerzon's velocity and energy
%   vectors, RV and RE, for sources placed at R and emitting signals G. For
%   K frequencies and N sources, G should be an K-by-N matrix and R should
%   be N-by-3.

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
%     [1] Gerzon (1992) General Metatheory of Auditory Localisation.
%     [2] Heller et al. (2008) Is My Decoder Ambisonic?

narginchk(2,2);

R = normalizeVector(R, 2); % normalize source directions
P = abs(G).^2; % source power

rV = real((G * R)./(sum(G,2) * ones(1,3)));
rE = (P * R)./(sum(P,2) * ones(1,3));

end