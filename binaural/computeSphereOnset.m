function D = computeSphereOnset(a,c,T,R,varargin)
%COMPUTESPHEREONSET Compute onset for a rigid-sphere HRIR.
%   D = COMPUTESPHEREONSET(A,C,T,R) computes the onset, D, of a rigid-
%   sphere HRIR of radius A (in m), given the speed of sound C (in m/s), 
%   angle of incidence T (in deg.), and non-dimensional source distance, R.
%   All inputs must be real scalars. The computed onset is returned in 
%   seconds. 
%
%   D = COMPUTESPHEREONSET(A,C,T,R,F) optionally includes a correction 
%   based on the specified maximum frequency F (in Hz). F is typically the
%   the chosen sampling frequency.

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
%   Copyright (c) 2020 Princeton University
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

% Refs:
%   [1]. Sridhar and Choueiri (2021) - Validation of the Minimum-Phase 
%   Representation of the Rigid-Sphere Head-Related Transfer Function.

narginchk(4,5)

if nargin < 5
    corrFlag = false;
else
    F = varargin{1};
    muStar = pi*F*a/c;
    cStar = c/(1+(0.5094*(2*muStar^2)^(-1/3)));
    corrFlag = true;
end

if R == inf
    if T <= 90
        D = -cosd(T)*(a/c);
    else
        if corrFlag
            D = deg2rad(T-90)*(a/cStar);
        else
            D = deg2rad(T-90)*(a/c);
        end
    end
else
    T0 = acosd(1/R);
    hatG = sqrt(R^2-1);
    G = sqrt(R^2-(2*R*cosd(T))+1);
    if T < T0
        D = G*(a/c);
    else
        if corrFlag
            D = hatG*(a/c) + deg2rad(T-T0)*(a/cStar);
        else
            D = (hatG + deg2rad(T-T0))*(a/c);
        end
    end
end

end
