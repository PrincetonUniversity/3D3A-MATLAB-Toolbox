function rP = wrap_wierstorf2013(b, Fs, varargin)
%WRAP_WIERSTORF2013 Wierstorf's azimuthal binaural localization model.
%   R = WRAP_WIERSTORF2013(B,FS,...) computes a predicted localization
%   vector (given in Cartesian coordinates) for binaural signals B at
%   sampling rate FS. The following additional inputs may also be
%   specified:
%
%   'model'     Localization model to use to estimate azimuth, either
%               'dietz2011' (default) or 'lindemann1986'.
%
%   'weights'   Frequency weighting to use, either 'no_spectral_weighting',
%               'rms_weighting' (default), or 'raatgever_weighting'.
%
%   'hrtf'      Path to .sofa file containing the listener's HRTFs.
%
%   'lookup'    ITD lookup table previously generated for the listener's
%               HRTFs.
%
%   NOTE: one of either 'hrtf' or 'lookup' must be specified. If the lookup
%   table is not provided, the HRTFs will be used to generate a new one.
%
%   See also WIERSTORF2013_GETLOOKUPTABLE, WIERSTORF2013_ESTIMATEAZIMUTH.

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2018 Princeton University
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

indx = find(strcmpi(varargin,'model'),1,'first');
if ~isempty(indx)
    modelToUse = varargin{indx+1};
else
    modelToUse = 'dietz2011';
end

indx = find(strcmpi(varargin,'weights'),1,'first');
if ~isempty(indx)
    specWeights = varargin{indx+1};
else
    specWeights = 'rms_weighting';
end

indx = find(strcmpi(varargin,'hrtf'),1,'first');
if ~isempty(indx)
    hrtfFile = varargin{indx+1};
else
    hrtfFile = '';
end

indx = find(strcmpi(varargin,'lookup'),1,'first');
if ~isempty(indx)
    % Import lookup table
    lookupTable = varargin{indx+1};
elseif ~isempty(hrtfFile)
    % Construct new lookup table from HRTF
    lookupTable = wierstorf2013_getLookupTable(hrtfFile, modelToUse);
else
    error('Must either provide lookup table or specify HRTF.');
end

if exist('wierstorf2013_estimateazimuth','file') == 2
    [AZ,~,~,~,~] = wierstorf2013_estimateazimuth(b, lookupTable, 'fs',...
        Fs, modelToUse, specWeights, 'remove_outlier');
else
    error('wierstorf2013_estimateazimuth from AMTOOLBOX not found.');
end

[rP(:,1),rP(:,2),rP(:,3)] = sph2cart(shiftdim(AZ)*pi/180,0,1);

end