function mag = getFiltMag(in,Fs,F)
%GETFILTMAG Compute magnitude of filter at specific frequency.
%   B = GETFILTMAG(A,FS) computes the magnitude, in dB, of the filter
%   with impulse response (IR), A, at the frequency at or closest to 1 kHz 
%   assuming a sampling rate of FS (in Hz).
%       If A is a vector, B is a scalar.
%       If A is a matrix, the individual IRs must be stored as columns and
%       the resulting B is a row vector.
%
%   B = GETFILTMAG(...,F) computes the magnitude at a frequency at or 
%   closest to F Hz.

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

narginchk(2,3);

if nargin < 3
    F = 1000;
end

in = shiftdim(in);
[inLen,~] = size(in);

fVec = getFreqVec(Fs,inLen);
[~,fIndx] = min(abs(fVec-F));
inTF = getMagSpecdB(in);
mag = inTF(fIndx,:);

end
