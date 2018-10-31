function [rI, D] = merimaa2005(B, varargin)
%MERIMAA2005 Merimaa's intensity vector and diffuseness parameter.
%   [I, D] = MERIMAA2005(X) computes the intensity vector and diffuseness
%   parameter using ambisonics signals X. The signals should be in columns
%   in ACN order. I will be a K-by-3 matrix, where K is the number of rows
%   in X, and D will be a K-by-1 vector.

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

%   References:
%     [1] Merimaa and Pulkki (2005) Spatial Impulse Response Rendering I: 
%         Analysis and Synthesis.

narginchk(1, inf);

rho_0 = 1.225; % density of air in kg/m^3
Z0 = rho_0 * getSoundSpeed(); % acoustic impedance of air in kg/(m^2 s)

% renormalize ambisonics signals if necessary
indx = find(strcmpi(varargin,'ambNorm'),1,'first');
if ~isempty(indx)
    B = ambReNormalize(B, varargin{indx+1}, 'N3D');
end

W = B(:,1)/sqrt(2);
Y = B(:,2)/sqrt(3);
Z = B(:,3)/sqrt(3);
X = B(:,4)/sqrt(3);

ex = [1 0 0];
ey = [0 1 0];
ez = [0 0 1];

Xp = X*ex + Y*ey + Z*ez;

temp = real(diag(conj(W))*Xp);
rI = (sqrt(2)/Z0)*temp;

% normalize the vector(s)
indx = find(strcmpi(varargin,'normalize'),1,'first');
if ~isempty(indx)
    rI = normalizeVector(rI,2);
end

F = sqrt(dot(temp,temp,2));
E = abs(W).^2 + dot(Xp,Xp,2)/2;
D = 1 - sqrt(2)*F./E;

end