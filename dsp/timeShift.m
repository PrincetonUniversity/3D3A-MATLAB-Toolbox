function shiftSignalMat = timeShift(inputSignalMat,shiftVec,METHOD)
%TIMESHIFT Delay/advance signals in time.
%   shiftSignalMat = TIMESHIFT(inputSignalMat,shiftVec) shifts signals
%   stored as the columns of inputSignalMat by the corresponding sample
%   amounts in the row vector, shiftVec. If shiftVec is a scalar, the same 
%   shift is applied to all columns of inputSignalMat. shiftVec can include 
%   fractional-sample values. The length of shiftVec should be the number 
%   of columns of inputSignalMat. shiftSignalMat has the same dimensions as
%   inputSignalMat.
%
%   shiftSignalMat = TIMESHIFT(...,METHOD) additionally specifies the shift
%   METHOD which can be linear ('lin') or circular ('circ'). If 'lin' is 
%   specified, the length of a shifted signal is the length of the
%   corresponding input signal plus the maximum absolute value of the shift,
%   rounded up to the next highest integer, specified in shiftVec. The
%   default is 'circ'.

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

narginchk(2,3);

if nargin < 3
    METHOD = 'circ';
end

inputSignalMat = shiftdim(inputSignalMat);
[nRows,nCols] = size(inputSignalMat);
if isscalar(shiftVec)
    shiftVec = shiftVec*ones(1,nCols);
end
shiftVec = reshape(shiftVec,[1 length(shiftVec)]);

switch lower(METHOD)
    case 'circ'
        phaseVec = exp(-1i*2*pi*getFreqVec(1,nRows)*shiftVec);
        shiftSignalMat = ifft(fft(inputSignalMat).*phaseVec,[],1,...
            'symmetric');
    case 'lin'
        maxShift = max(ceil(abs(shiftVec)));
        paddedLen = 2^nextpow2(nRows+maxShift);
        padSigMat = padarray(inputSignalMat,[paddedLen-nRows 0],0,'post');
        phaseVec = exp(-1i*2*pi*getFreqVec(1,paddedLen)*shiftVec);
        shiftSignalMat = ifft(fft(padSigMat).*phaseVec,[],1,...
            'symmetric');
        shiftSignalMat = shiftSignalMat(1:nRows+maxShift,:);
    otherwise
        error('Invalid method specification.')
end

end

