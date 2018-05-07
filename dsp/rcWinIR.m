function outputIR = rcWinIR(inputIR,winLen,r,winStart,oLenOpts)
%RCWINIR Apply raised-cosine window to impulse response(s).
%   outputIR = RCWINIR(inputIR,winLen) applies a raised-cosine window  of 
%   length winLen starting at the first sample of each impulse response in 
%   inputIR. If inputIR contains more than one IR, the IRs must be stored
%   as columns. The raised-cosine window onset and roll-off are both set to
%   0.5 (corresponding to half the winLen). For more on the raised-cosine
%   window parameters, see RAISEDCOSWIN.
%
%   outputIR = RCWINIR(...,r) specifies r, which is a 1x2 vector where the 
%   first element is the fraction of winLen for the window onset and the 
%   second element is the fraction of winLen for the window roll-off. 
%   The default value for r is [0.5,0.5].
%
%   outputIR = RCWINIR(...,winStart) specifies the starting sample of the
%   window. If inputIR contains more than one IR, winStart can be a vector
%   whose length is equal to the number of columns in inputIR, each
%   containing a starting sample value for the corresponding IR. The
%   default value for winStart is 1.
%
%   outputIR = RCWINIR(...,oLenOpts) specifies two options for the length
%   of outputIR. The two options are: 'full' - IRs in outputIR have same
%   length as IRs in inputIR, and 'trunc' - IRs in outputIR have
%   length winLen (default).
%
%   See also RAISEDCOSWIN, RECTWINIR, TUKEYWINIR.

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
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

narginchk(2,5);

inputIR = shiftdim(inputIR);
[irLen,numCols] = size(inputIR);

if ~isscalar(winLen)
    error('winLen must be a scalar.')
end

if nargin < 5
    oLenOpts = 'trunc';
end

if nargin < 4
    winStart = 1;
end

if nargin < 3
   r = [0.5,0.5]; 
end

if isscalar(winStart)
    winStart = winStart*ones(numCols,1);
else
    winStart = shiftdim(winStart);
    if length(winStart) ~= numCols
        error('Length of winStart must equal number of columns in inputIR')
    end
end

switch lower(oLenOpts)
    case 'full'
        outputIR = zeros(irLen,numCols);
    case 'trunc'
        outputIR = zeros(winLen,numCols);
    otherwise
        error('Invalid specification for oLenOpts')
end

winVec = raisedCosWin(winLen,r);
for ii = 1:numCols
    outputIR(1:winLen,ii) = inputIR(winStart(ii):(winStart(ii)+winLen-1),...
        ii).*winVec;
end

end

