function S = normalizeVector(R, DIM)
%NORMALIZEVECTOR Normalize vectors to unit magnitude.
%   S = NORMALIZEVECTOR(R) returns a unit vector S that points in the same
%   direction as the vector R.
%
%   S = NORMALIZEVECTOR(M) normalizes the vectors stored in the columns of
%   the matrix M.
%
%   S = NORMALIZEVECTOR(M, DIM) normalizes the vectors in matrix M along
%   dimension DIM.

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
%   Copyright (c) 2017 Princeton University
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
%   =======================================================================

if isvector(R)
    S = R/norm(R);
elseif ismatrix(R)
    if nargin < 2
        DIM = 1;
    elseif DIM > 2
        error('Invalid dimension.');
    end
    
    switch DIM
        case 1
            S = R*diag(1./sqrt(dot(R,R,DIM)));
        case 2
            S = diag(1./sqrt(dot(R,R,DIM)))*R;
    end
else
    error('Must be matrix or vector.');
end

end