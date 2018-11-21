function [E, fc] = scharer2009(h, h0, Fs, fc)
%SCHARER2009 Scharer's spectral errors in auditory bands.
%   C = SCHARER2009(X, [], FS, FC) computes errors in the energy spectrum
%   of X (i.e., abs(fft(X)).^2) in auditory bands with center frequencies
%   specified by FC. Auditory bands are modeled by gammatone filters, as
%   implemented in the LTFAT toolbox.
%   
%   C = SCHARER2009(X, X0, FS, FC) computes errors relative to the energy
%   spectrum of X0.
%
%   [C, FC] = SCHARER2009(...) returns the center frequencies also.
%   
%   Suggested center frequencies for the equivalent rectangular bandwidth
%   (ERB) scale:
%       fc = erbspacebw(f_low,f_high); % Needs LTFAT toolbox.
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
%     [1] Scharer and Lindau (2009) Evaluation of Equalization Methods for
%         Binaural Signals.

FFTLen = size(h,1);
H = abs(fft(h,FFTLen,1));
if isempty(h0)
    H0 = ones(size(h));
else
    if ~isequal(size(h),size(h0))
        error('Reference and input signals must be the same length.');
    end
    H0 = abs(fft(h0,FFTLen,1));
end

g = getGammatoneFilters(fc, Fs, FFTLen);
Gmag = abs(fft(g,FFTLen,1)); % FFTLen-by-length(fc)

E = mag2db(sqrt(((Gmag.')*(abs(H).^2))./((Gmag.')*(abs(H0).^2)))); % Eq. (9)

end