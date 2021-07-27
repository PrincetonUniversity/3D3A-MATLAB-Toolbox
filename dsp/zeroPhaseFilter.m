function y = zeroPhaseFilter(varargin)
%ZEROPHASEFILTER Zero-phase filtering. 
%   Y = ZEROPHASEFILTER(H,X) filters the signal(s) in X with the FIR filter 
%   H to create filtered signal(s) Y. If X is a matrix, each column is 
%   treated as a separate signal. The dimensions of Y will match those of 
%   X. The signal(s) in X are assumed to be causal.
% 
%   Y = ZEROPHASEFILTER(SOS,G,X) filters the signal(s) in X with the 
%   second-order section (SOS) filter described by the matrix SOS and the
%   vector G.
%
%   Needs: Signal Processing Toolbox.
%
%   See also FILTFILT, FFTCONV.

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
%   Copyright (c) 2021 Princeton University
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

narginchk(2,3);

switch nargin
    case 2
        convertFlag = false;
        h = varargin{1};
        x = varargin{2};
    case 3
        convertFlag = true;
        sos = varargin{1};
        g = varargin{2};
        x = varargin{3};
end

flipFlag = isrow(x);
if flipFlag
    x = x.';
end

xLen = size(x,1);

if convertFlag
    if ~bitget(abs(xLen),1) % xLen is even
        diracVec = [1;zeros(xLen,1)];
        hLen = xLen+1;
    else % xLen is odd
        diracVec = [1;zeros(xLen-1,1)];
        hLen = xLen;
    end
    
    h = g*sosfilt(sos,diracVec);
    H = getMagSpec(h);
else
    h = shiftdim(h);
    hLen = length(h);
    H = getMagSpec(h);
    if ~bitget(abs(hLen),1) % hLen is even
        H = [H;abs(H(hLen))];
        hLen = hLen+1;
    end
end

h_lp = circshift(ifft(H,'symmetric'),(hLen-1)/2);
filtx = circshift(fftConv(h_lp,x,'lin'),-(hLen-1)/2);
y = filtx(1:xLen,:);

if flipFlag
    y = y.';
end

end
