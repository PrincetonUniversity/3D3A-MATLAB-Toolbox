function onsetMat = thresholdIRs(irData,thp,varargin)
%THRESHOLDIRS Find the onset(s) of input IR(s) using thresholding.
%   D = THRESHOLDIRS(X) identifies the onset samples D for input IRs X
%   using a 10% threshold. If X is an matrix M-by-N, thresholding is
%   performed along its columns, and D is then a 1-by-N row vector. If X is
%   an M-by-N-by-K 3D array, thresholding is performed along the third
%   dimension and D is then an N-by-M matrix.
%
%   Similarly, if X is a 2D cell array of size M-by-N, then D will be an
%   M-by-N matrix. Note that each element of the cell array must be a
%   vector, and thresholding will be performed along each of those vectors.
%
%   D = THRESHOLDIRS(X,THP) specifies the threshold value to use to compute
%   the onsets. THP takes on values between 0 and 1 and corresponds to a
%   fraction of the maximum absolute value of each IR. The first sample in
%   a given IR to exceed THP times the maximum absolute amplitude of that
%   IR is considered the onset sample of that IR.
%
%   D = THRESHOLDIRS(X,THP,Name1,Value1,...) specifies optional
%   comma-separated pairs of Name,Value arguments, where Name is the 
%   argument name and Value is the corresponding value. Name must appear 
%   inside single quotes (' '). You can specify the following:
%
%   'resample'      Resample X prior to computing onsets. This value must
%                   be a scalar > 0, where a value in the range (0,1)
%                   corresponds to downsampling, and a value > 1
%                   corresponds to upsampling. Onsets are returned at the
%                   original sampling rate.
%
%   Needs: Signal Processing Toolbox.
%
%   See also ESTIMATEIRONSET.

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

% Needs at least 1 input argument
narginchk(1,4);

% Use 10% (-20 dB) threshold by default
if nargin == 1
    thp = 0.1;
end

% Check input data type
array2D = (~iscell(irData) && ismatrix(irData));
array3D = (~iscell(irData) && ndims(irData)==3);
cell2D = (iscell(irData) && ismatrix(irData));
if ~(array2D || array3D || cell2D)
    error(['Input data must be either: a 2D numeric array, a 3D',...
        ' numeric array, or a 2D cell array.']);
end
    
% Get array dimensions
switch ndims(irData)
    case 2
        [nRows,nCols] = size(irData);
    case 3
        [nRows,nCols,~] = size(irData);
end

% Perform resampling ('upsample' included for backwards-compatibility)
indx = find(strcmpi(varargin,'resample')|strcmpi(varargin,'upsample'),1);
if ~isempty(indx)
    [p,q] = rat(varargin{indx+1});
    if array2D
        irData = resample(irData,p,q); % From Signal Processing Toolbox
    elseif array3D
        tempIR = resample(irData(1,1,:),p,q);
        resampleData = zeros(nRows,nCols,length(tempIR));
        for ii = 1:nRows
            for jj = 1:nCols
                resampleData(ii,jj,:) = resample(irData(ii,jj,:),p,q);
            end
        end
        irData = resampleData;
    elseif cell2D
        for ii = 1:nRows
            for jj = 1:nCols
                tempIR = irData{ii,jj};
                if ~isvector(tempIR)
                    error('Cell array elements must be vectors.')
                end
                irData{ii,jj} = resample(shiftdim(tempIR),p,q);
            end
        end
    end
else
    p = 1;
    q = 1;
end

% Perform thresholding
if array2D
    onsetMat = zeros(1,nCols);
    for jj=1:nCols
        peakVal = max(abs(irData(:,jj)));
        onsetMat(jj) = find(abs(irData(:,jj)) >= (thp*peakVal),1,...
            'first')*(q/p);
    end
elseif array3D
    onsetMat = zeros(nRows,nCols);
    for ii=1:nRows
        for jj=1:nCols
            peakVal = max(abs(irData(ii,jj,:)));
            onsetMat(ii,jj) = find(abs(irData(ii,jj,:)) >= ...
                (thp*peakVal),1,'first')*(q/p);
        end
    end
elseif cell2D
    onsetMat = zeros(nRows,nCols);
    for ii=1:nRows
        for jj=1:nCols
            peakVal = max(abs(irData{ii,jj}));
            onsetMat(ii,jj) = find(abs(irData{ii,jj}) >= (thp*peakVal),...
                1,'first')*(q/p);
        end
    end
end

end
