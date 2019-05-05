function [SPL,fVec] = getISO226Values(p,f)
%GETISO226VALUES Equal-loudness contours as described in ISO226
%   [SPL,F] = GETISO226VALUES(PHON) returns dB SPL values at 29 frequencies
%   (as specified in IS0226) that defines the equal-loudness contour for a 
%   specified PHON value (in dB), where 1 PHON = 1 dB SPL at 1 kHz. Also
%   returned is F, the vector of frequencies at which the contour is
%   defined. The frequencies span from 20 Hz to 12.5 kHz and are returned
%   as a column vector in units of Hz. The SPL values of the equal-loudness
%   contour are also returned as a column vector. PHON must be a real
%   scalar between 0 and 90.
%
%   [SPL,F] = GETISO226VALUES(PHON,FV) returns the equal-loudness contours
%   sampled at frequencies specified in the vector FV. The frequencies in
%   FV must be specified in Hz. The SPL values of the equal-loudness
%   contours are computed using cubic spline interpolation. The returned
%   vector of frequencies, F, corresponds to the vector of frequencies at
%   which the contour is originally specified.

%   This function is based on a function written by Jeff Tackett on 
%   03/01/05, the original copyright notice for which is provided below:

%   =======================================================================
%   ORIGINAL COPYRIGHT NOTICE
%
%   Copyright (c) 2009, Jeff Tackett
%   All rights reserved.
% 
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions are
%   met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the 
%       distribution
% 
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
%   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
%   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
%   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
%   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
%   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
%   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
%   (INCLUDING NEGLIGENCE OR OTHERWISE)ARISING IN ANY WAY OUT OF THE USE OF 
%   THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%   =======================================================================

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

narginchk(1,2);

validateattributes(p,{'double'},{'scalar','nonempty','nonnan',...
    'finite','real','nonnegative','<=',90},'getISO226Values','PHON',1)

% Values tabulated in ISO226
fVec = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 ...
     1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500].';

af = [0.532 0.506 0.480 0.455 0.432 0.409 0.387 0.367 0.349 0.330 0.315 ...
      0.301 0.288 0.276 0.267 0.259 0.253 0.250 0.246 0.244 0.243 0.243 ...
      0.243 0.242 0.242 0.245 0.254 0.271 0.301].';

Lu = [-31.6 -27.2 -23.0 -19.1 -15.9 -13.0 -10.3 -8.1 -6.2 -4.5 -3.1 ...
       -2.0  -1.1  -0.4   0.0   0.3   0.5   0.0 -2.7 -4.1 -1.0  1.7 ...
        2.5   1.2  -2.1  -7.1 -11.2 -10.7  -3.1].';

Tf = [ 78.5  68.7  59.5  51.1  44.0  37.5  31.5  26.5  22.1  17.9  14.4 ...
       11.4   8.6   6.2   4.4   3.0   2.2   2.4   3.5   1.7  -1.3  -4.2 ...
       -6.0  -5.4  -1.5   6.0  12.6  13.9  12.3].';

% Setup user-defined values for equation
Ln = p;

% Deriving sound pressure level from loudness level (See ISO226 Sec. 4.1)
Af = 4.47E-3 * (10.^(0.025*Ln) - 1.15) + (0.4*10.^(((Tf+Lu)/10)-9 )).^af;
SPL = ((10./af).*log10(Af)) - Lu + 94;

if nargin > 1
    validateattributes(f,{'double'},{'vector','nonempty','nonnan',...
        'finite','real','nonnegative'},'getISO226Values','FV',2)
    SPL = spline(fVec,SPL,f);
end

end
