function az = vec2azim(r)
%VEC2AZIM Compute the azimuthal direction of a vector.
%   AZ = VEC2AZIM(R) computes the azimuthal direction AZ (in radians) of
%   the vector R. If R is an M-by-3 matrix, AZ will be a length M column
%   vector. If R is a cell array, AZ will be a cell array of the same size;
%   each element of R must be a matrix, and each corresponding element of
%   AZ will be a vector.
%
%   See also CART2SPH, CART2POL, AZIM2VEC.

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolkit.
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

% if r is a cell array, recursively loop over each element
if iscell(r)
    az = cell(size(r));
    for ii = 1:numel(r)
        az{ii} = azim2vec(r{ii});
    end
    return
end

% r must be a matrix
if ~ismatrix(r)
    error('Input must be a matrix (or a cell array of matrices)!');
end

switch size(r,2)
    case 2
        [az,~] = cart2pol(r(:,1),r(:,2));
    case 3
        [az,~,~] = cart2sph(r(:,1),r(:,2),r(:,3));
    otherwise
        error('Vector dimensions must be either 2 or 3.')
end

end