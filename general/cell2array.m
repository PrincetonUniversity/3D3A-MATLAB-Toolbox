function M = cell2array(C,OPT)
%CELL2ARRAY Convert a nested cell array to a multi-dimensional array.
%   M = CELL2ARRAY(C) converts the nested cell array, C, to a
%   multi-dimensional array, M. This command may also be specified as:
%       M = CELL2ARRAY(C,'no-reduce');
%
%   M = CELL2ARRAY(C,OPT) specifies whether or not singleton dimensions 
%   should be removed from the multi-dimensional array. Specifying OPT as 
%   'reduce' results in the removal of any singleton dimensions in M.
%   Otherwise, singleton dimensions are not removed in the conversion.

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

narginchk(1,2);

if nargin < 2
    OPT = 'no-reduce';
end

[rowDims,colDims] = getCellDimensions(C);
dimVec = reshape([rowDims;colDims],1,[]);
dimVecLen = length(dimVec);
tempMat2 = [];
for ii = 1:rowDims(1)
    tempMat1 = [];
    for jj = 1:colDims(1)
        if iscell(C{ii,jj})
            tempMat1 = cat(dimVecLen-1,tempMat1,cell2array(C{ii,jj}));
        else
            tempMat1 = cat(dimVecLen-1,tempMat1,C{ii,jj});
        end
    end
    tempMat2 = cat(dimVecLen,tempMat2,tempMat1);
    M = permute(tempMat2,[dimVecLen,dimVecLen-1,1:(dimVecLen-2)]);
end

if strcmpi(OPT,'reduce')
    M = squeeze(M);
end

end
