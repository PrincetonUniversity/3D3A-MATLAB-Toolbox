function d = computeDirectionalError(r1,r2,errorType)
%COMPUTEDIRECTIONALERROR Directional error between two vectors.
%   E = COMPUTEDIRECTIONALERROR(R1,R2) computes the angular error E between
%   vectors (given in Cartesian coordinates) R1 and R2.
%
%   E = COMPUTEDIRECTIONALERROR(R1,R2,TYPE) computes error of the following
%   types:
%       'l2'    l2 norm (Euclidean distance) between R1 and R2 after
%               normalizing each vector (also 'error').
%
%       'angle' Angular distance between R1 and R2 (also 'angular').

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
    errorType = 'angle';
end

if iscell(r1) && iscell(r2)
    % call this function recursively for cell array inputs
    d = cell(size(r1));
    for ii = 1:numel(r1)
        d{ii} = computeDirectionalError(r1{ii},r2{ii},errorType);
    end
else
    r1n = normalizeVector(r1,2);
    r2n = normalizeVector(r2,2);
    
    switch lower(errorType)
        case {'l2','error'}
            temp = r1n - r2n;
            d = sqrt(dot(temp,temp,2));
        case {'angle','angular'}
            d = real(acos(dot(r1n,r2n,2)));
    end
    d(isnan(d)) = 0;
end

end