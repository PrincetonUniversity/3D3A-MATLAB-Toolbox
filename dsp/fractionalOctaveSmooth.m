function Hsm = fractionalOctaveSmooth(H, varargin)
%FRACTIONALOCTAVESMOOTH Fractional-octave smoothing of transfer functions.
%   HS = FRACTIONALOCTAVESMOOTH(H,N) computes the 1/N-octave smoothed
%   transfer function HS given the raw transfer function H. If H is a
%   matrix, the smoothing operation is carried out along its columns.
%
%   HS = FRACTIONALOCTAVESMOOTH(H,M) uses a precomputed smoothing matrix M.
%
%   HS = FRACTIONALOCTAVESMOOTH(H,N,METHOD) specifies the smoothing METHOD
%   to be used. By default, 'tylka' is assumed.
%
%   HS = FRACTIONALOCTAVESMOOTH(H,N,METHOD,WINTYPE) specifies the type of
%   the smoothing window, either 'rectangular' (default) or 'hanning'.
%
%   HS = FRACTIONALOCTAVESMOOTH(H,N,METHOD,WINTYPE,SCALE) performs
%   smoothing on the specified SCALE, either 'power' (default), 'dB',
%   'complex', or 'equiv-complex'.
%
%   See also COMPUTESMOOTHINGMATRIX.

%   ==============================================================================
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
%     [1] Hatziantoniou and Mourjopoulos (2000) Generalized
%         Fractional-Octave Smoothing of Audio and Acoustic Responses.
%     [2] Tylka et al. (2017) A Generalized Method for Fractional-Octave
%         Smoothing of Transfer Functions that Preserves Log-Frequency
%         Symmetry.

narginchk(2,5);

FFTLen = size(H,1);

if ~isscalar(varargin{1})
    % Use precomputed smoothing matrix
    M = varargin{1};
    if FFTLen ~= size(M,2)
        error('Size mismatch between input transfer function and smoothing matrix.');
    end
else
    % Compute smoothing matrix using specified parameters
    FRAC = varargin{1};
    PARAMS = {FRAC, 'tylka','rectangular','power'};
    for ii = 2:numel(varargin)
        PARAMS{ii} = varargin{ii};
    end
    METHOD = PARAMS{2};
    WINTYPE = PARAMS{3};
    SCALE = PARAMS{4};
    M = computeSmoothingMatrix(FFTLen, FRAC, METHOD, WINTYPE);
end

switch(lower(SCALE))
    case 'power'
        Hsm1 = sqrt(M*(abs(H).^2));
    case 'db'
        Hsm1 = db2mag(M*mag2db(abs(H)));
    case 'complex'
        Hsm1Re = M*real(H);
        Hsm1Im = M*imag(H);
        Hsm1 = Hsm1Re+1i*Hsm1Im;
    case 'equiv-complex'
        Hsm1Re = M*real(H);
        Hsm1Im = M*imag(H);
        Hsm1a = Hsm1Re+1i*Hsm1Im;
        Hsm1b = sqrt(M*(abs(H).^2));
        Hsm1 = (Hsm1a.*Hsm1b)./(abs(Hsm1a));
end


Hsm = [Hsm1; flip(conj(Hsm1(2:end-1,:)))];

end