function h = getPattersonFilters(fc,Fs,IRLen,cFlag)
%GETPATTERSONFILTERS Auditory filters that represent human critical bands.
%   H = GETPATTERSONFILTERS(FC,FS,N) or H = GETPATTERSONFILTERS(FC,FS,N,...
%   'zerophase') computes length N, zero-phase impulse responses (IRs) of 
%   Patterson's auditory-band filters with center frequencies specified by
%   the vector FC and given at a sampling rate FS. The output, H, will be 
%   an N-by-M matrix of IRs, where M is the number of elements in FC.
%
%   H = GETPATTERSONFILTERS(...,'causal') makes the IRs causal by applying 
%   a shift equal to half the length of the IRs.
%
%   See also GETGAMMATONEFILTERS.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2019 Princeton University
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
%     [1] Salomons (1995) Coloration and Binaural Decoloration of Sound due
%         to Reflections.

narginchk(3,4);

if nargin < 4
    cFlag = 'zerophase';
end

if ~isrow(fc)
    fc = shiftdim(fc).';
end
numfc = length(fc);
f = getFreqVec(Fs, IRLen);
specLen = ceil((1+IRLen)/2);

Wfc = fc2erb(fc,2); % Salomons, Eq. 5.11
temp = 4*abs(f(1:specLen)*ones(1,numfc) - ...
    ones(specLen,1)*fc)./(ones(specLen,1)*Wfc);
H = (1 + temp).*exp(-temp); % Salomons, Eq. 5.9

h_zp = ifft(H,IRLen,1,'symmetric');

switch lower(cFlag)
    case 'zerophase'
        h = h_zp;
    case 'causal'
        h = shiftSignal(h_zp,IRLen/2);
    otherwise
        error('Fourth input may be ''zerophase'' or ''causal'' only.')
end

end
