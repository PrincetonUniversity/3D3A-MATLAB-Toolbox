function [outputData,outputPos] = sortData(inputData,posMat,DIM,ROUND)
%SORTDATA Sort inputData based on sorting of posMat.
%   [outputData,outputPos] = SORTDATA(inputData,posMat) removes any
%   non-singleton dimensions in posMat, sorts the first column in ascending
%   order, and correspondingly sorts inputData. inputData and posMat must 
%   have the same number of rows. The sorted data is returned and has the
%   same dimensions as inputData. The sorted posMat is also returned.
%
%   [outputData,outputPos] = SORTDATA(...,DIM) additionally specifies 
%   specific columns to sort sequentially. For example, DIM = [2,1] first 
%   sorts by column 2, then by column 1.
%
%   [outputData,outputPos] = SORTDATA(...,ROUND) additionally specifies 
%   the number of digits to round data in posMat prior to sorting. For more
%   information, see ROUND.
%
%   See also SORTROWS, ROUND.

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
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
%   ==============================================================================

narginchk(2,4);

if nargin < 4
    ROUND = 3;
end

if nargin < 3
    DIM = 1;
end

[outputPos,sortIndices] = sortrows(round(posMat,ROUND),DIM);
outputData = inputData(:,sortIndices);

end