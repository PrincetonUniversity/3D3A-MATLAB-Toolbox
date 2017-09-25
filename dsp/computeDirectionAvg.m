function avgVal = computeDirectionAvg(inputSpectra,METHOD)
%COMPUTEDIRECTIONAVG Average input spectra over "directions" (columns).
%   avgVal = COMPUTEDIRECTIONAVG(inputSpectra) computes the linear average,
%   over "directions" (columns), of the spectra in inputSpectra. 
%   inputSpectra may be a vector or 2D matrix. If it's a matrix, individual
%   spectra must be stored as columns and avgVal is a column vector with 
%   length equal to the number of rows in inputSpectra. If inputSpectra is
%   a vector, then it is forced to a row vector and avgVal is a scalar.
%
%   avgVal = COMPUTEDIRECTIONAVG(...,METHOD) additionally specifies the 
%   method used to average the data in inputSpectra across columns to 
%   produce avgVal. The options are 'lin' (linear average, default), and 
%   'rms' (root-mean-square).
%
%   See also MEAN, RMS, COMPUTESPECTRUMAVG.

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

narginchk(1,2);

if nargin < 2
    METHOD = 'lin';
end

if isvector(inputSpectra)
    inputSpectra = reshape(inputSpectra,[1,length(inputSpectra)]);
end

switch lower(METHOD)
    case 'lin'
        avgVal = mean(inputSpectra,2);
    case 'rms'
        avgVal = rms(inputSpectra,2);
    otherwise
        error('Invalid input for METHOD');
end

end

