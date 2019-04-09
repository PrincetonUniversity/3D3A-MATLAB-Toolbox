function [alignedIRs,dMat] = timeAlignIRs(irData,thp)
% TIMEALIGNIRS Time align input impulse responses
%   [Y,D] = TIMEALIGNIRS(X) identifies the onset samples in X using
%   thresholding with a relative threshold of 10% of the absolute maximum
%   and uses this to align the IRs. If X is a matrix, it is assumed the
%   columns are individual IRs. The matrix of aligned IRs, Y, and the
%   shift, D, in samples, used to align each IR to the IR with the smallest
%   onset are returned.
%
%   [Y,D] = TIMEALIGNIRS(X,THP) optionally specifies the threshold value to
%   use. THP can take values between 0 and 1 with 1 corresponding to a
%   relative threshold of 100% of the absolute maximum.
%
%   See also THRESHOLDIRS.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   Joseph G. Tylka <josephgt@princeton.edu>
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

narginchk(1,2);

if nargin < 2
    thp = 0.1; % Use 10% (-20 dB) threshold by default
end

switch ndims(irData)
    case 2
        [nRows,nCols] = size(irData);
    case 3
        [nRows,nCols,IRLen] = size(irData);
    otherwise
        error('X has incompatible input dimensions.')
end
onsetMat = thresholdIRs(irData,thp);
dMat = onsetMat-min(min(onsetMat));

if iscell(irData)
    alignedIRs = cell(nRows,nCols);
    for ii = 1:nRows
        for jj = 1:nCols
            alignedIRs{ii,jj} = circshift(irData{ii,jj},[-dMat(ii,jj) 0]);
        end
    end
else
    switch ndims(irData)
        case 2
            alignedIRs = zeros(nRows,nCols);
            for ii = 1:nCols
                alignedIRs(:,ii) = circshift(irData(:,ii),[-dMat(1,ii) 0]);
            end
        case 3
            alignedIRs = zeros(nRows,nCols,IRLen);
            for ii = 1:nRows
                for jj = 1:nCols
                    alignedIRs(ii,jj,:) = circshift(irData(ii,jj,:),...
                        [0 0 -dMat(ii,jj)]);
                end
            end
        otherwise
            error('X has incompatible input dimensions.')
    end
end

end
