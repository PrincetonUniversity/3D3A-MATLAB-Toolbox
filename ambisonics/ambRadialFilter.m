function d = ambRadialFilter(l,k,r)
%AMBRADIALFILTER Radial distance coding filters for near-field ambisonics.
%   D = AMBRADIALFILTER(L,K,R) computes the ambisonics radial distance
%   coding filter D, for angular wavenumber K and order L, for a source
%   distance R.
%
%   L, K, and R may be matrices, in which case they must all have the same
%   size or else one or more must be a scalar. In this case, D will have
%   the same dimensions as the matrices.
%
%   Alternatively, K may be a column vector and R may be a row vector, in
%   which case L must be a scalar or match K or R, and D will be LENGTH(K)-
%   by-LENGTH(R).
%
%   See also AMBPOINTSOURCE.

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

[l,k,kr] = prepInputs(l,k,r);
d = ((1i.^(l+1)).*k).*sphericalHankelH(l,1,kr);
d(kr==0) = +~l;

end

function [l,k,kr] = prepInputs(l,k,r)

leqk = all(size(l)==size(k));
leqr = all(size(l)==size(r));
keqr = all(size(k)==size(r));
scalars = [isscalar(l) isscalar(k) isscalar(r)];

if sum(scalars) >= 2 % at least two of l, k, and r are scalars
    kr = k.*r;
    return
end % else, at most one scalar

if isscalar(l) && keqr % l is the only scalar and k and r are the same size
    kr = k.*r;
    return
end

if isscalar(k) && leqr % k is the only scalar and l and r are the same size
    kr = k*r;
    return
end

if isscalar(r) && leqk % r is the only scalar and l and k are the same size
    kr = k*r;
    return
end

if leqk && leqr && keqr % l, k, and r are all the same size
    kr = k.*r;
    return
end

if iscolumn(k) && isrow(r) % k*r is a nonscalar matrix
    kr = k*r;
    k = repmat(k,1,length(r)); % k has same size as kr
    
    % now try to make l the same size as kr
    if isscalar(l)
        l = repmat(l,size(kr));
        return
    elseif leqk % l is a column vector like k
        l = repmat(l,1,length(r));
        return
    elseif leqr % l is a row vector like r
        l = repmat(l,length(k),1);
        return
    end
end

error(['Could not combine input matrices. '...
    'Two or more of L, K, and R must be scalars or matrices of the same size. '...
    'Alternatively, K may be a column-vector and R a row-vector, '...
    'in which case L must be a scalar or match K or R.']);
    
end