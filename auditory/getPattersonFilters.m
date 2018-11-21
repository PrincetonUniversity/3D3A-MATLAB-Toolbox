function h = getPattersonFilters(fc, Fs, IRLen)
%GETPATTERSONFILTERS Auditory filters that represent human critical bands.
%   H = GETPATTERSONFILTERS(FC,FS,N) computes length N impulse responses of
%   Patterson's auditory-band filters with center frequencies specified by
%   the vector FC and given at a sampling rate FS.
%   
%   The output, H, will be an N-by-M matrix of impulse responses, where M
%   is the number of elements in FC.
%
%   See also GETGAMMATONEFILTERS.

%   =======================================================================
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

%   References:
%     [1] Salomons (1995) Coloration and Binaural Decoloration of Sound due
%         to Reflections.

numfc = length(fc);
f = getFreqVec(Fs, IRLen);
specLen = 1 + IRLen/2;

Wfc = fc2erb(fc,2); % Salomons, Eq. 5.11
temp = 4*abs(f(1:specLen)*ones(1,numfc) - ones(specLen,1)*fc)./(ones(specLen,1)*Wfc);
H = (1 + temp).*exp(-temp); % Salomons, Eq. 5.9

h = ifft(H,IRLen,1,'symmetric');

end