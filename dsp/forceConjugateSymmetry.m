function Hsym = forceConjugateSymmetry(H)
%FORCECONJUGATESYMMETRY Ensure a conjugate-symmetric transfer function.
%   HS = FORCECONJUGATESYMMETRY(H) returns HS, the conjugate-symmetric
%   version of the transfer function H. H may be a row or column vector, or
%   a matrix, where each column is a different transfer function. This
%   function may be used instead of specifying the 'symmetric' option when
%   computing an IFFT.
%
%   If H is a cell array, FORCECONJUGATESYMMETRY is performed recursively
%   on each element of H.
%
%   See also IFFT.

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

if iscell(H)
    Hsym = cell(size(H));
    for ii = 1:numel(H)
        Hsym{ii} = forceConjugateSymmetry(H{ii});
    end
end

rowVec = isrow(H);
if rowVec
    H = H.';
end

HLen = size(H,1);

Hsym = H;
Hsym(1,:) = real(H(1,:)); % DC is real
if ~mod(HLen,2) % HLen is even
    Hsym(HLen/2+1,:) = real(H(HLen/2+1,:)); % Nyquist is real
    Hsym(HLen/2+2:HLen,:) = conj(flipud(H(2:HLen/2,:)));
else % Hlen is odd
    Hsym((HLen+1)/2+1:HLen,:) = conj(flipud(H(2:(HLen+1)/2,:)));
end

if rowVec
    Hsym = Hsym.';
end

end
