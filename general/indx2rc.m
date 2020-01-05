function varargout = indx2rc(X,DIM)
%INDX2RC Row and columns indices from linear index.
%   [R,C] = INDX2RC(X,DIM) takes a vector of indices X and determines the
%   row and column number assuming the indices correspond to positions in a
%   matrix, M, of dimensions given by DIM. DIM must be a 2-element vector
%   specifying the number of rows and columns of M. The returned values, R
%   and C, correspond to the row and column number, respectively, 
%   corresponding to each index in X. R and C will both be column vectors
%   with the same length as X.
%
%   Y = INDX2RC(X,DIM) returns a 2-column matrix Y with as many rows as 
%   elements in X. The first column of Y is R and the second column is C.
%
%   Example:
%       A = magic(3);
%       A([2,4,9]) = -1;
%       X = find(A == -1);
%       [R,C] = INDX2RC(X,[3,3]);
%       Y = INDX2RC(X,[3,3]);

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2019 Princeton University
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

% Check input count
narginchk(2,2);

% Validate required inputs
validateattributes(X,{'numeric'},{'vector','nonempty','nonnan',...
    'finite','integer'},'indx2rc','X',1)
validateattributes(DIM,{'numeric'},{'vector','nonempty','nonnan',...
    'finite','integer','numel',2},'indx2rc','DIM',2)

X = shiftdim(X);
DIM = shiftdim(DIM);

C = ceil(X/DIM(2));
R = mod(X,DIM(1));
R(R == 0) = DIM(1);

switch nargout
    case 1
        varargout = {[R,C]};
    case 2
        varargout = {R,C};
    otherwise
        error('Only 1 or 2 outputs may be requested.')
end

end
