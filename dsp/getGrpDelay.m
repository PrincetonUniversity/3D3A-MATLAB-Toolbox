function [grpDelaySpec,avgGrpDelay] = getGrpDelay(inputIR,fS,AVGRANGE)
%GETGRPDELAY Compute the group delay of a transfer function given its 
%impulse response.
%   [grpDelaySpec,avgGrpDelay] = GETGRPDELAY(inputIR,fS) computes the group 
%   delay of the transfer function with impulse response given in inputIR. 
%   The outputs are specified in samples. inputIR may be a vector or 
%   matrix. If inputIR is a matrix, the IRs must be stored as columns and 
%   the output avgGrpDelay will be a row vector while grpDelaySpec will 
%   have the same dimensions as inputIR. fS is the sampling rate in Hz.
%
%   [grpDelaySpec,avgGrpDelay] = GETGRPDELAY(...,AVGRANGE) optionally
%   specifies the frequency range over which averaging of grpDelaySpec 
%   should be performed to compute avgGrpDelay. AVGRANGE must be specified 
%   as a row vector [w1,w2] containing frequencies specified in Hz such 
%   that w1 < w2.

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
%
%   Reference:
%       [1]. https://ccrma.stanford.edu/~jos/filters/Numerical_Computation_
%            Group_Delay.html#14646

narginchk(2,3);

inputIR = shiftdim(inputIR);
FFTLen = size(inputIR,1);

if nargin < 3
    AVGRANGE = [0,fS/2];
end

freqVec = getFreqVec(fS,FFTLen);
[~,lIndx] = min(abs(freqVec-AVGRANGE(1)));
[~,hIndx] = min(abs(freqVec-AVGRANGE(2)));

inputTF = fft(inputIR);

% Replace very small values in inputTF by 1 prior to division.
minmag = 10*eps;
skipIndxs = (abs(inputTF) < minmag);
inputTF(skipIndxs) = 1;

% See [1] for theory behind following calculation.
grpDelaySpec = real(fft(diag(0:FFTLen-1)*inputIR)./inputTF);

% Replace samples where calculation would have "blown up" by 0.
grpDelaySpec(skipIndxs) = 0;

% Compute average group delay.
avgGrpDelay = mean(grpDelaySpec(lIndx:hIndx,:),1,'omitnan');

end
