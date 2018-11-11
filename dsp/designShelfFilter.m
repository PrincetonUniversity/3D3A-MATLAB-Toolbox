function [b, a] = designShelfFilter(G, Wc, Q, type)
%DESIGNSHELFFILTER Low- and high-shelf filters.
%   [B,A] = DESIGNSHELFFILTER(G,WC,Q,TYPE) returns filter coefficients B
%   and A for a shelf filter of type TYPE with normalized corner frequency
%   WC = FC/FS/2, shelf gain G dB, and quality factor Q.

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

K = tan(pi * Wc / 2);
V0 = db2mag(G);
Z = 1/Q;

if(V0 < 1)
    V0 = 1/V0;
end

if(( G > 0 ) && (strcmp(type,'low')))
    b0 = (1 + sqrt(V0)*Z*K + V0*K^2) / (1 + Z*K + K^2);
    b1 =             (2 * (V0*K^2 - 1) ) / (1 + Z*K + K^2);
    b2 = (1 - sqrt(V0)*Z*K + V0*K^2) / (1 + Z*K + K^2);
    a1 =                (2 * (K^2 - 1) ) / (1 + Z*K + K^2);
    a2 =             (1 - Z*K + K^2) / (1 + Z*K + K^2);
elseif (( G < 0 ) && (strcmp(type,'low')))
    b0 =             (1 + Z*K + K^2) / (1 + Z*sqrt(V0)*K + V0*K^2);
    b1 =                (2 * (K^2 - 1) ) / (1 + Z*sqrt(V0)*K + V0*K^2);
    b2 =             (1 - Z*K + K^2) / (1 + Z*sqrt(V0)*K + V0*K^2);
    a1 =             (2 * (V0*K^2 - 1) ) / (1 + Z*sqrt(V0)*K + V0*K^2);
    a2 = (1 - Z*sqrt(V0)*K + V0*K^2) / (1 + Z*sqrt(V0)*K + V0*K^2);
elseif (( G > 0 ) && (strcmp(type,'high')))
    b0 = (V0 + Z*sqrt(V0)*K + K^2) / (1 + Z*K + K^2);
    b1 =             (2 * (K^2 - V0) ) / (1 + Z*K + K^2);
    b2 = (V0 - Z*sqrt(V0)*K + K^2) / (1 + Z*K + K^2);
    a1 =              (2 * (K^2 - 1) ) / (1 + Z*K + K^2);
    a2 =           (1 - Z*K + K^2) / (1 + Z*K + K^2);
elseif (( G < 0 ) && (strcmp(type,'high')))
    b0 =               (1 + Z*K + K^2) / (V0 + Z*sqrt(V0)*K + K^2);
    b1 =                  (2 * (K^2 - 1) ) / (V0 + Z*sqrt(V0)*K + K^2);
    b2 =               (1 - Z*K + K^2) / (V0 + Z*sqrt(V0)*K + K^2);
    a1 =             (2 * ((K^2)/V0 - 1) ) / (1 + Z/sqrt(V0)*K + (K^2)/V0);
    a2 = (1 - Z/sqrt(V0)*K + (K^2)/V0) / (1 + Z/sqrt(V0)*K + (K^2)/V0);
else
    b0 = V0;
    b1 = 0;
    b2 = 0;
    a1 = 0;
    a2 = 0;
end

a = [  1, a1, a2];
b = [ b0, b1, b2];

end