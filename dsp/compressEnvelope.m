function b = compressEnvelope(a,C)
%COMPRESSENVELOPE Compress the envelope of an input signal.
%   B = COMPRESSENVELOPE(A) compresses the envelope of an input signal, A, 
%   using the algorithm described by Bernstein et al. [1]. The compression 
%   applied is typically used to simulate peripheral auditory compression 
%   by the basilar membrane. If A is a matrix of signals, each column of A
%   is treated as a separate signal and compressed independently.
%
%   B = COMPRESSENVELOPE(A,C) optionally specifies, C, the compression 
%   factor to use in the compression algorithm. C can take values between 0 
%   and 1. The default value is 0.2.
%
%   Needs: Signal Processing Toolbox.

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
%       [1]. Bernstein et al. (1999) - The normalized interaural 
%       correlation: Accounting for NoS? thresholds obtained with Gaussian 
%       and ?low-noise? masking noise.

narginchk(1,2);

if nargin < 2
    C = 0.2;
end

% Check inputs
validateattributes(a,{'double'},{'2d','nonempty','nonnan','finite'},...
    'applyPeripheralCompression','A',1);
validateattributes(C,{'double'},{'scalar','nonempty','nonnan','finite',...
    'real','nonnegative','<=',1},'applyPeripheralCompression','C',2);

% Compute (upper) Hilbert envelope
a_mean = mean(a);
a_hilb = hilbert(a_mean);
a_e = abs(a_hilb);

% Compute signal with compressed envelope
b = ((a_e.^(C-1)).*(a-a_mean)) + a_mean;

end
