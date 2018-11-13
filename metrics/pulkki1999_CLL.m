function [CLL, fc] = pulkki1999_CLL(b, Fs, fc)
%PULKKI1999_CLL Pulkki's composite loudness level model.
%   [CLL,FC] = PULKKI1999_CLL(B,FS) computes the composite loudness level
%   CLL at each center frequency FC given binaural signals B at sampling
%   rate FS.
%
%   [CLL,~] = PULKKI1999_CLL(B,FS,FC) uses specifed center frequencies FC.
%
%   See also GETGAMMATONEFILTERS, BUTTER, FFTCONV, FILTER.

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

%   References:
%     [1] Pulkki et al. (1999) Analyzing Virtual Sound Source Attributes
%         Using a Binaural Auditory Model.
%     [2] Huopaniemi et al. (1999) Objective and Subjective Evaluation of
%         Head-Related Transfer Function Filter Design.

IRLen = size(b,1);

if nargin < 3
    % Compute center frequencies of the Gammatone filters
    N = 42;
    f1 = 200;
    QBW = 9.26449 * 24.7;
    fc = logspace(log10(f1+QBW), log10(Fs/2+QBW), N+1) - QBW;
    fc = fc(1:end-1);
end

pinkNoiseSource = dsp.ColoredNoise(1, IRLen, 1);
x = shiftdim(step(pinkNoiseSource));
bxL = fftConv(b(:,1),x,'lin');
bxR = fftConv(b(:,2),x,'lin');

h = getGammatoneFilters(fc, Fs, IRLen);
[B, A] = butter(1,1000/(Fs/2));

% Apply Gammatone filters
xl = fftConv(real(h),bxL,'lin');
xr = fftConv(real(h),bxR,'lin');
% Rectify signal
xl = xl.*(sign(xl)+1)/2;
xr = xr.*(sign(xr)+1)/2;
% Apply low-pass filter
xl = filter(B,A,xl,[],1);
xr = filter(B,A,xr,[],1);
% Estimate loudness
LL = sqrt(sqrt(sum(xl .* xl,1)/size(xl,1)));
LR = sqrt(sqrt(sum(xr .* xr,1)/size(xr,1)));

CLL = 10*log2(LL + LR) + 40;

end