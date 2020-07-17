function [grpDelaySpec,avgGrpDelay] = getGrpDelay(inputIR,Fs,varargin)
%GETGRPDELAY Group delay from impulse response (IR).
%   [S,M] = GETGRPDELAY(A,Fs) computes the group delay of the IR(s) in A 
%   and returns the group delay spectrum, S, and the average group delay, 
%   M.
%       If A is a vector, S will be a column vector with the same length as 
%       A, and M will be a scalar.
%       If A is a matrix with dimensions P-by-N, S will have the same size 
%       as A while M will be a row vector of length N.
%   The group delay values in S and M are specified in samples. The
%   sampling rate, Fs, must be specified in Hz.
%
%   [S,M] = GETGRPDELAY(A,Fs,AVGRANGE) optionally specifies the frequency 
%   range over which averaging of S should be performed to compute M. 
%   AVGRANGE must be specified as a row vector [FL,FU] containing 
%   frequencies specified in Hz such that 0 <= FL <= FU <= Fs/2. The 
%   default value is [0,Fs/2].
%
%   See also ESTIMATEIRONSET.

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

%   Reference:
%       [1]. https://ccrma.stanford.edu/~jos/filters/Numerical_Computation_
%            Group_Delay.html#14646

% Check input count
narginchk(2,3);

% Verify required inputs
validateattributes(inputIR,{'numeric'},{'2d','nonempty','nonnan',...
    'finite'},'getGrpDelay','A',1)
validateattributes(Fs,{'numeric'},{'scalar','real','finite','nonnan',...
    'positive'},'getGrpDelay','Fs',2)

inputIR = shiftdim(inputIR); % If inputIR is a vector, force to a column.
irLen = size(inputIR,1);

% Verify optional inputs
if nargin < 3
    avgRange = [0,Fs/2];
else
    avgRange = varargin{1};
    validateattributes(avgRange,{'numeric'},{'vector','real',...
        'finite','nonnan','nonnegative','numel',2,'<=',Fs/2},...
        'getGrpDelay','AVGRANGE',3)
end

if avgRange(1) >= avgRange(2)
    error(['When AVGRANGE is specified as [FL,FU], FU must be greater',...
        ' than FL.'])
end

% === Begin main calculations ===

freqVec = getFreqVec(Fs,irLen);
[~,lIndx] = min(abs(freqVec-avgRange(1)));
[~,hIndx] = min(abs(freqVec-avgRange(2)));

inputTF = fft(inputIR);

% Replace very small values in inputTF by 1 prior to division.
minmag = 10*eps;
skipIndxs = (abs(inputTF) < minmag);
inputTF(skipIndxs) = 1;

% See [1] for theory behind following calculation.
grpDelaySpec = real(fft(diag(0:irLen-1)*inputIR)./inputTF);

% Replace samples where above calculation would have "blown up" by 0.
grpDelaySpec(skipIndxs) = 0;

% Compute average group delay.
avgGrpDelay = mean(grpDelaySpec(lIndx:hIndx,:),1,'omitnan');

end
