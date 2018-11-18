function d = ambRadialFilter(l,k,r)
%AMBRADIALFILTER Radial distance coding filters for near-field ambisonics.
%   D = AMBRADIALFILTER(L,K,R) computes the ambisonics radial distance
%   coding filter D, for angular wavenumber K and order L, for a source
%   distance R.
%
%   If K is a vector, D will be a LENGTH(K) vector in the same orientation
%   as K. R must be a scalar.

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

if ~isscalar(l)
    error('L must be a scalar.');
end

if ~isvector(k)
    error('K must be a vector.');
end

if ~isscalar(r)
    error('R must be a scalar.');
end

if k(1)==0
    DIM = find(size(k) ~= 1, 1, 'first');
    k = k(2:end);
    dropZero = true;
    zeroVal = +~l;
else
    dropZero = false;
end

d = ((1i^(l+1)).*k).*sphericalHankelH(l,1,k*r);

if dropZero
    if ~isempty(DIM)
        d = cat(DIM,zeroVal,d);
    elseif dropZero && isempty(DIM)
        d = zeroVal;
    end
end

end