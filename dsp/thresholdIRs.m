function onsetMat = thresholdIRs(irData,thp,varargin)
%THRESHOLDIRs Find the onset(s) of input IR(s) using thresholding.
%   onsetMat = THRESHOLDIRs(irData) identifies the onset samples in input 
%   IRs using a 10% threshold. If irData is a matrix, the input IRs must be
%   stored as columns. onsetMat is then a row vector. If irData is a
%   tensor (3D matrix), the deepest dimension is assumed to contain the IR
%   data. onsetMat is then a matrix. For example, if irData has dimensions 
%   2-by-3-by-100, it is assumed that 100 is the length of each of the IRs.
%
%   onsetMat = THRESHOLDIRs(...,thp) specifies the threshold value to use
%   to compute the onset. thp takes on values between 0 and 1 and
%   corresponds to a fraction of the maximum absolute value of each IR. The
%   first sample in a given IR to exceed thp times the maximum absolute
%   amplitude of that IR is considered the onset sample of that IR.
%
%   onsetMat = THRESHOLDIRs(...,Name1,Value1,...) specifies optional 
%   comma-separated pairs of Name,Value arguments, where Name is the 
%   argument name and Value is the corresponding value. Name must appear 
%   inside single quotes (' '). You can specify the following:
%
%   'resample'      Resample irData prior to computing onsets. This must
%                   be a scalar > 0 such that a value in the range (0,1)
%                   corresponds to downsampling, and a value > 1
%                   corresponds to upsampling. This corresponds to
%                   specifying parameter 'P' in the 'resample' function
%                   assuming 'Q' is set to 1.

%   ==============================================================================
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

% Needs at least 1 input argument
if nargin == 0
    onsetMat = 0;
    return
end

% Use 10% (-20 dB) threshold by default
if nargin == 1
    thp = 0.1;
end

% Perform resampling ('upsample' included for backwards-compatibility)
indx = findInCell(varargin,'resample')+findInCell(varargin,'upsample');
if indx && ~iscell(irData) && ~(ndims(irData)-2)
    p = varargin{indx+1};
    irData = resample(irData,p,1);
else
    p = 1;
end

switch ndims(irData)
    case 2
        [nRows,nCols] = size(irData);
    case 3
        [nRows,nCols,~] = size(irData);
end

if iscell(irData)
    onsetMat = zeros(nRows,nCols);
    for ii=1:nRows
        for jj=1:nCols
            peakVal = max(abs(irData{ii,jj}));
            onsetMat(ii,jj) = find(abs(irData{ii,jj}) >= (thp*peakVal),1,...
                'first');
        end
    end
else
    switch ndims(irData)
        case 2
            onsetMat = zeros(1,nCols);
            for jj=1:nCols
                peakVal = max(abs(irData(:,jj)));
                onsetMat(jj) = find(abs(irData(:,jj)) >= (thp*peakVal),1,...
                    'first')/p;
            end
        case 3
            onsetMat = zeros(nRows,nCols);
            for ii=1:nRows
                for jj=1:nCols
                    peakVal = max(abs(irData(ii,jj,:)));
                    onsetMat(ii,jj) = find(abs(irData(ii,jj,:)) >= ...
                        (thp*peakVal),1,'first');
                end
            end
    end
end

end