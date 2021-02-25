function [y,lag] = xcoh(x1,x2,varargin)
%XCOH Cross-coherence
%   Y = XCOH(X1,X2) computes the cross-coherence between input signals X1 
%   and X2 and returns the output in Y. Cross-coherence is a normalized
%   cross-correlation between X1 and X2 and measures the similarity between
%   X1 and shifted (i.e., lagged) copies of X2.
%       If X1 and X2 are vectors, they must have the same length, N. Y will 
%       be a column vector with length N.
%       If X1 and X2 are matrices, they must have the same dimensions and
%       the cross-coherence of corresponding columns is evaluated.
%
%   Y = XCOH(X1,X2,METHOD) optionally allows the method used to compute
%   cross-coherence to be specified. METHOD can take the following two
%   options:
%       1. 'weighted' (default) - Cross-coherence is computed as the
%       cross-correlation of X1 and X2 normalized by the geometric mean of
%       the total energies of X1 and X2. This corresponds to a weighted sum
%       of phase errors.
%
%       2. 'unweighted' - Cross-coherence is computed as the 
%       cross-correlation of X1 and X2 deconvolved by the convolution of 
%       zero-phase versions of X1 and X2. This corresponds to an unweighted 
%       sum of phase errors.
%
%       3. 'timedomain' - Cross-coherence is computed in the time domain.
%       The output should be very close to, if not exactly the same as,
%       that obtained when using the 'weighted' option.
%
%   Y = XCOH(X1,X2,METHOD,FRANGE) optionally specifies the frequency range
%   over which the cross-coherence is computed. FRANGE must be specified as
%   a two-element vector, [F1,F2], where F1 and F2 are the lower- and
%   upper-bound frequencies. F1 and F2 must be specified such that 0
%   corresponds to DC and 1 to the Nyquist frequency. This option is not
%   applicable when METHOD is set to 'timedomain'.
%
%   [Y,L] = XCOH(...) additionally returns lag indices, L.
%
%   Needs: Signal Processing Toolbox.
%
%   See also XCORR.

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
%   Copyright (c) 2020 Princeton University
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

%   Ref:
%       [1]. Nam et al. (2008) - On the Minimum-Phase Nature of Head-
%       Related Transfer Functions.

% Check input count
narginchk(2,5);

% Validate required inputs
validateattributes(x1,{'double'},{'2d','finite','nonnan','nonempty'},...
    'xcoh','X1',1)
x1 = shiftdim(x1); % If vector, force to column
validateattributes(x2,{'double'},{'2d','finite','nonnan','nonempty'},...
    'xcoh','X2',2)
x2 = shiftdim(x2); % If vector, force to column

if size(x1) ~= size(x2)
    error(['X1 and X2 must have the same length (if vectors) or',...
        ' dimensions (if matrices).'])
end

% Validate optional inputs
if nargin < 5 % Hidden input
    matchXCorrFlag = false;
else
    matchXCorrFlag = varargin{3};
    validateattributes(matchXCorrFlag,{'logical'},{'scalar','nonempty'},...
        'xcoh','matchXCorrFlag (hidden option)',5)
end

if nargin < 4
    frange = [0,1];
else
    frange = varargin{2};
    validateattributes(frange,{'numeric'},{'vector','nonempty','real',...
        'nonnegative','<=',1},'xcoh','FRANGE',4)
end

if nargin < 3
    method = 'weighted';
else
    method = varargin{1};
    validateattributes(method,{'char'},{'scalartext','nonempty'},...
        'xcoh','METHOD',3)
end

numRows = size(x1,1);
switch lower(method)
    % 'weighta' and 'weightb' options included for backwards compatibility
    case {'weighted','unweighted','weighta','weightb'}
        if matchXCorrFlag
            mxl = numRows-1;
            ceilLog2 = nextpow2(2*numRows-1);
            m2 = 2^ceilLog2;
            X1 = fft(x1,m2);
            X2 = fft(x2,m2);
            numRows = size(X1,1);
        else
            X1 = fft(x1);
            X2 = fft(x2);
        end

        % Compute normalized frequency indices
        fVec = linspace(0,2-(2/numRows),numRows);
        [~,fL] = min(abs(fVec-frange(1)));
        [~,fU] = min(abs(fVec-frange(2)));
        X_mask = zeros(size(X1));
        X_mask(fL:fU,:) = 1;
        X_mask = forceConjugateSymmetry(X_mask);
        X1_mask = X1.*X_mask;
        X2_mask = X2.*X_mask;
    case 'timedomain'
        numCols = size(x1,2);
        y = zeros(2*numRows-1,numCols);
    otherwise
        error('Invalid METHOD specification.')
end

switch lower(method)
    % 'weighta' option included for backwards compatibility
    case {'weighted','weighta'}
        if isreal(x1) && isreal(x2)
            y_un = ifft(X1_mask.*conj(X2_mask),'symmetric');
        else
            y_un = ifft(X1_mask.*conj(X2_mask));
        end
        x1_sq = sqrt(sum(abs(X1_mask).^2).*sum(abs(X2_mask).^2))*...
            (1/numRows);
        if matchXCorrFlag
            y1 = y_un./repmat(x1_sq,numRows,1);
            y = [y1(m2 - mxl + (1:mxl),:); y1(1:mxl+1,:)];
            lag = -mxl:mxl;
        else
            y = y_un./repmat(x1_sq,numRows,1);
            lag = 0:(numRows-1);
        end
    % 'weightb' option included for backwards compatibility
    case {'unweighted','weightb'}
        if isreal(x1) && isreal(x2)      
            y_un = ifft((abs(X_mask).^2).*exp(1i*angle(X1.*conj(X2))),...
                'symmetric');
        else
            y_un = ifft((abs(X_mask).^2).*exp(1i*angle(X1.*conj(X2))));
        end
        x1_sq = sum(abs(X_mask).^2)*(1/numRows);
        if matchXCorrFlag
            y1 = y_un./repmat(x1_sq,numRows,1);
            y = [y1(m2 - mxl + (1:mxl),:); y1(1:mxl+1,:)];
            lag = -mxl:mxl;
        else
            y = y_un./repmat(x1_sq,numRows,1);
            lag = 0:(numRows-1);
        end
    case 'timedomain'
        for ii = 1:numCols
            [y(:,ii),lag] = xcorr(x1(:,ii),x2(:,ii),'coeff');
        end
    otherwise
        error('Invalid METHOD specification.')
end

end
