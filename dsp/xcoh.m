function [y,lag] = xcoh(x1,x2)
%XCOH Cross-coherence
%   Y = XCOH(X1,X2) computes the cross-coherence between input signals X1 
%   and X2 and returns the output in Y.
%       If X1 and X2 are vectors, they must have the same length ,N. Y will 
%       be a column vector with length (2*N) - 1.
%       If X1 and X2 are matrices, they must have the same dimensions and
%       the cross-coherence of corresponding columns is evaluated.
%
%   [Y,L] = XCOH(...) additionally returns lag indices, L.
%
%   Needs: Signal Processing Toolbox.
%
%   See also XCORR.

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
%       [1]. Nam et al. (2008) - On the Minimum-Phase Nature of Head-
%       Related Transfer Functions.

narginchk(2,2);

validateattributes(x1,{'double'},{'2d','finite','nonnan','nonempty'},...
    'xcoh','X1',1)
x1 = shiftdim(x1);
validateattributes(x2,{'double'},{'2d','finite','nonnan','nonempty'},...
    'xcoh','X2',2)
x2 = shiftdim(x2);

if size(x1) ~= size(x2)
    error(['X1 and X2 must have the same length (if vectors) or',...
        ' dimensions (if matrices).'])
end

[numRows,numCols] = size(x1);
lag = zeros(numCols,(2*numRows)-1);
y = zeros((2*numRows)-1,numCols);
for ii = 1:numCols
    [y_xcorr,lag(ii,:)] = xcorr(x1(:,ii),x2(:,ii));
    % Take absolute value below to allow complex X1 and X2.
    y(:,ii) = y_xcorr/sqrt(sum(abs(x1(:,ii)).^2)*sum(abs(x2(:,ii)).^2));
end

end
