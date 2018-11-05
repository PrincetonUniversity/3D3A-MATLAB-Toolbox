function ncc = computeNCC(ref,test)
%COMPUTENCC Compute maximum of the normalized cross-correlation.
%   ncc = COMPUTENCC(ref,test) computes the maximum of the normalized
%   cross-correlation (ncc) between ref and test. Both ref and test may be 
%   vectors, or matrices with the same dimensions. If matrices, it is 
%   assumed that the signals are stored as columns. ncc ranges from -1 to 
%   1.

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

%   References:
%       [1]. Xie (2013) - Head-Related Transfer Function and Virtual
%       Auditory Display, Second Edition, J. Ross. 

narginchk(2,2);

validateattributes(ref,{'double'},{'2d','nonempty','nonnan','finite',...
    'ndims',2},'computeNCC','ref',1);
validateattributes(test,{'double'},{'2d','nonempty','nonnan','finite',...
    'size',size(ref)},'computeNCC','test',2);

ref = shiftdim(ref); % If row vector, convert to column vector.
test = shiftdim(test);

% The following formulas were adapted from page 138 of Xie [1].
numCols = size(ref,2);
ncc = zeros(1,numCols);
for ii = 1:numCols
    currentXCorr = xcorr(ref(:,ii),test(:,ii),'coeff');
    [~,nccIndx] = max(abs(currentXCorr));
    ncc(ii) = currentXCorr(nccIndx);
end

end
