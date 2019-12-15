function [y,l] = getInterpMax(x,DIM)
%GETINTERPMAX Quadratically-interpolated maximum
%   Y = GETINTERPMAX(X) returns the quadratically-interpolated maximum of
%   the discrete-time sequence, X, following the algorithm specified by J.
%   O. Smith [1].
%       If X is a vector, GETINTERPMAX operates along the non-singleton 
%       dimension and returns Y as a scalar.
%       If X is a matrix, GETINTERPMAX operates along columns and returns Y
%       as a row vector.
%
%   Y = GETINTERPMAX(X,DIM) specifies the dimension along which the maximum
%   should be computed.
%
%   [Y,L] = GETINTERPMAX(...) also returns, in L, the index/indices 
%   corresponding to the maximum values.
%
%   See also MAX.

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

%   Ref:
%       [1]. J.O. Smith (2011) - Spectral Audio Signal Processing, 
%       http://ccrma.stanford.edu/~jos/sasp/

narginchk(1,2)

if nargin < 2
    [beta,beta_indx] = max(x,[],'omitnan');
else
    [beta,beta_indx] = max(x,[],DIM,'omitnan');
end

x = shiftdim(x);
numCols = length(beta);
y = zeros(1,numCols);
l = zeros(1,numCols);
xLen = size(x,1);
for ii = 1:numCols
    if beta_indx(ii) == 1 || beta_indx(ii) == xLen
        l(ii) = beta_indx(ii);
        y(ii) = beta(ii);
    else
        gamma = x(beta_indx(ii)+1,ii);
        alpha = x(beta_indx(ii)-1,ii);
        p = 0.5*((alpha-gamma)/(alpha-(2*beta(ii))+gamma));
        l(ii) = beta_indx(ii) + p;
        y(ii) = beta(ii)-(0.25*(alpha-gamma)*p);
    end 
end

end
