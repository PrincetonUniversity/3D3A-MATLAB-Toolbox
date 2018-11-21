function r = azim2vec(az, n)
%AZIM2VEC Convert an azimuthal direction to a unit vector.
%   R = AZIM2VEC(AZ) computes a unit vector R pointing in the direction AZ
%   (given in radians). If AZ is a vector, R will be an M-by-3 matrix,
%   where M = LENGTH(AZ). If AZ is a cell array, R will be a cell array of
%   the same size; each element of AZ must be a vector, and each
%   corresponding element of R will be a matrix.
%
%   R = AZIM2VEC(AZ,N) specifies the dimension of each row of R, either 2
%   or 3 (default).
%
%   See also SPH2CART, POL2CART, VEC2AZIM.

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

% by default, create 3D vectors
if nargin < 2
    n = 3;
end

% if az is a cell array, recursively loop over each element
if iscell(az)
    r = cell(size(az));
    for ii = 1:numel(az)
        r{ii} = azim2vec(az{ii});
    end
    return
end

% az must be a vector
if ~isvector(az)
    error('Input must be a vector (or a cell array of vectors)!');
end

r = zeros(numel(az),n);

switch n
    case 2
        [r(:,1), r(:,2)] = pol2cart(az(:),1);
    case 3
        [r(:,1), r(:,2), r(:,3)] = sph2cart(az(:),0,1);
    otherwise
        error('Vector dimensions must be either 2 or 3.')
end

end