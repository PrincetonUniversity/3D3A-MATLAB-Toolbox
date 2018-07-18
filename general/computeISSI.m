function ISSI = computeISSI(inputTFMat,ISSIType)
%COMPUTEISSI Compute the inter-sweet-spot isolation between a pair of
%sweet-spots.
%   ISSI = COMPUTEISSI(inputTFMat,ISSIType) computes the inter-sweet-spot
%   isolation between pairs of sweet-spots using transfer function data for
%   each pair of those sweet-spots. inputTFMat must be either a N-by-4 or
%   N-by-8 array depending on the ISSIType specification as described
%   below, where N corresponds to the number of sample points in a given
%   transfer function.
%
%   =============================================
%   ISSIType                inputTFMat Dimensions
%   =============================================
%   'Driver Out'            N-by-4
%   'Driver In'             N-by-4
%   'Passenger Out'         N-by-4
%   'Passenger In'          N-by-4
%   'Combined Out'          N-by-8
%   'Combined In'           N-by-8
%
%   The output, ISSI, is a vector of length N specified in dB.

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

%   References:
%       [1] Gallian (2017) - Independent Audio Program Delivery to Multiple
%       Listeners in Reverberating Environments.pdf

narginchk(2,2);

[~,numTFs] = size(inputTFMat);
if (strcmpi(ISSIType,'driver out') || strcmpi(ISSIType,'driver in') ...
        || strcmpi(ISSIType,'passenger out') || ...
        strcmpi(ISSIType,'passenger in')) && numTFs ~= 4
    error('inputTFMat dimensions do not correspond to specified ISSIType.')
elseif (strcmpi(ISSIType,'combined out') || ...
        strcmpi(ISSIType,'combined in')) && numTFs ~= 8
    error('inputTFMat dimensions do not correspond to specified ISSIType.')
end

switch lower(ISSIType)
    case 'driver out'
        H11 = inputTFMat(:,1);
        H12 = inputTFMat(:,2);
        H13 = inputTFMat(:,3);
        H14 = inputTFMat(:,4);
        ISSI = 10*log10((abs(H11).^2 + abs(H12).^2)./...
            (abs(H13).^2 + abs(H14).^2));
    case 'driver in'
        H21 = inputTFMat(:,1);
        H22 = inputTFMat(:,2);
        H23 = inputTFMat(:,3);
        H24 = inputTFMat(:,4);
        ISSI = 10*log10((abs(H21).^2 + abs(H22).^2)./...
            (abs(H23).^2 + abs(H24).^2));
    case 'passenger out'
        H41 = inputTFMat(:,1);
        H42 = inputTFMat(:,2);
        H43 = inputTFMat(:,3);
        H44 = inputTFMat(:,4);
        ISSI = 10*log10((abs(H43).^2 + abs(H44).^2)./...
            (abs(H41).^2 + abs(H42).^2));
    case 'passenger in'
        H31 = inputTFMat(:,1);
        H32 = inputTFMat(:,2);
        H33 = inputTFMat(:,3);
        H34 = inputTFMat(:,4);
        ISSI = 10*log10((abs(H33).^2 + abs(H34).^2)./...
            (abs(H31).^2 + abs(H32).^2));
    case 'combined out'
        ISSIDOut = computeISSI(inputTFMat(:,1:4),'driver out');
        ISSIPOut = computeISSI(inputTFMat(:,5:8),'passenger out');
        ISSI = min([ISSIDOut'; ISSIPOut'])';
    case 'combined in'
        ISSIDIn = computeISSI(inputTFMat(:,1:4),'driver in');
        ISSIPIn = computeISSI(inputTFMat(:,5:8),'passenger in');
        ISSI = min([ISSIDIn'; ISSIPIn'])';
    otherwise
        error('Invalid ISSIType specification.')
end

end

