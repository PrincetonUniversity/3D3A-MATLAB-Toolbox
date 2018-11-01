function h = getGammatoneFilters(fc, Fs, IRLen)
%GETGAMMATONEFILTERS Auditory filters that represent human critical bands.
%   h = GETGAMMATONEFILTERS(fc,Fs,N) computes length N impulse responses of
%   gammatone filters with center frequencies specified by the vector fc.
%   
%   This function is essentially a wrapper for the gammatonefir function in
%   the LTFAT toolbox. The output, h, will be an N-by-M matrix of filters,
%   where M is the number of elements in fc, rather than a length M cell
%   array.
%   
%   Needs LTFAT toolbox.
%
%   See also GAMMATONEFIR.

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

numfc = length(fc);
h = zeros(IRLen, numfc);

hcell = gammatonefir(fc, Fs, IRLen);

for ii = 1:numfc
    tempLen = length(hcell{ii}.h);
    h(1:tempLen,ii) = hcell{ii}.h;
end

% % Example plot of gammatone filters
%
% Fs = 96000;
% N = 4096;
% fc = erbspacebw(200,16000);
% h = getGammatoneFilters(fc, Fs, N);
% f = getFreqVec(Fs, N);
% Hmag = mag2db(abs(fft(h,N,1)));
% figure()
% hold all
% for ii = 1:length(fc)
% plot(f, Hmag(:,ii));
% end
% set(gca,'XScale','log')
% xlim([fc(1)/1.5, fc(end)*1.5])
% hold off

end