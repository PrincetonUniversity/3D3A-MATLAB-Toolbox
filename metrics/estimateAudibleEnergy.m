function E = estimateAudibleEnergy(x, Fs, FRANGE)
%ESTIMATEAUDIBLEENERGY Average energy in auditory critical bands.
%   E = ESTIMATEAUDIBLEENERGY(X,FS) computes the total audible energy E for
%   a signal X given at sample rate FS. The total audible energy is
%   computed as the sum of the signal energies in critical bands, modeled
%   here by ERB-spaced (equivalent rectangular bandwidth) gammatone filters.
%
%   E = ESTIMATEAUDIBLEENERGY(X,FS,[FL,FH]) computes the total audible
%   energy in the frequency range [FL,FH]. The default is [20,20000].
%
%   See also GETGAMMATONEFILTERS.

%   ==============================================================================
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

if nargin < 3
    FRANGE = [20,20000];
end

if exist('erbspacebw','file') == 2
    fc = erbspacebw(min(FRANGE),max(FRANGE));
else
    warning('erbspacebw from LTFAT not found, using 1/3-octave bands instead.');
    numfc = round(log2(max(FRANGE)/min(FRANGE))/(1/3));
    fc = logspace(log10(min(FRANGE)),log10(min(FRANGE)*2^(numfc/3)),numfc);
end

FFTLen = size(x,1);

h = getGammatoneFilters(fc, Fs, FFTLen);
Hmag = getMagSpec(h,1); % FFTLen-by-length(fc)

X = fft(x,FFTLen,1);
E = mean(((Hmag.')*(abs(X).^2)));

end