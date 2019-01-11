function out = applyFiltersInSeries(filtMat,in)
%APPLYFILTERSINSERIES Filter input signal by multiple filters sequentially.
%   C = APPLYFILTERSINSERIES(A,B) takes an input signal, B, and filters it 
%   by each FIR filter in A, sequentially, to produce the output, C. The
%   individual filters in A must be stored as columns. The input B can be a
%   matrix of multiple signals, each stored as a column. In this case, each
%   signal in B is filtered by the filters in A. The output, C, has the
%   same dimensions as the input, B.

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

narginchk(2,2);

% Check inputs
validateattributes(filtMat,{'double'},{'2d','nonempty','nonnan',...
    'finite'},'applyFiltersInSeries','A',1);
validateattributes(in,{'double'},{'2d','nonempty','nonnan',...
    'finite'},'applyFiltersInSeries','B',2);

numFilters = size(filtMat,2);
numSignals = size(in,2);
out = zeros(size(in)); % Initialize output matrix
for ii = 1:numSignals % Loop over each input signal
    currentSignal = in(:,ii);
    for jj = 1:numFilters % Apply filters sequentially (i.e. in series)
        currentSignal = filter(filtMat(:,jj),1,currentSignal);
    end
    out(:,ii) = currentSignal;
end

end
