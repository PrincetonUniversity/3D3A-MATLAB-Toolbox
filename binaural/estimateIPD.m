function varargout = estimateIPD(hL,hR,Fs,varargin)
%ESTIMATEIPD Estimate the interaural phase delay (IPD) spectrum in seconds.
%   IPD = ESTIMATEIPD(HL,HR,FS) estimates IPD for input HRIRs HL and HR at 
%   sampling rate FS specified in Hz. Negative IPDs correspond to sound 
%   sources on the left. 
%       If HL and HR are vectors, they must have the same dimensions and 
%       IPD will be a vector with these dimensions.
%       If HL and HR are matrices, the individual IRs must be stored as 
%       columns and the matrices must have the same dimensions. The
%       returned IPD will be a matrix with these dimensions.
%
%   [IPD,ITD] = ESTIMATEIPD(HL,HR,FS) also returns a frequency-independent
%   ITD in seconds by averaging the IPD spectrum over the entire frequency 
%   range.
%       If HL and HR are vectors, ITD will be a scalar.
%       If HL and HR are matrices, ITD will be a row vector with as many
%       elements as the number of columns in HL and HR.
%
%   [IPD,ITD] = ESTIMATEIPD(HL,HR,FS,[FL,FU]) computes ITD by averaging IPD
%   over the frequency range [FL,FU], where both frequencies are specified
%   in Hz. FL and FU must both be non-negative and less than the Nyquist
%   frequency corresponding to FS. Also, FL must be less than or equal to
%   FU.
%
%   [IPD,ITD] = ESTIMATEIPD(...,'linearfit') applies a linear fit to the
%   IPD spectrum before averaging to compute ITD.
%
%   See also ESTIMATEITD.

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

% === CHECK AND PARSE INPUTS; INITIALIZE REQUIRED VARIABLES === %

narginchk(3,5);
nargoutchk(1,2);

validateattributes(hL,{'numeric'},{'2d','nonempty','nonnan','finite'},...
    'estimateIPD','hL',1);
validateattributes(hR,{'numeric'},{'2d','nonempty','nonnan','finite',...
    'size',size(hL)},'estimateIPD','hR',2);
validateattributes(Fs,{'numeric'},{'scalar','nonempty','nonnan',...
    'finite','real','positive'},'estimateIPD','Fs',3);

if isrow(hL)
    % Convert row vectors to column vectors
    hL = shiftdim(hL); 
    hR = shiftdim(hR);
    tpFlag = true; % Flag to determine whether IPD should be transposed
else
    tpFlag = false;
end

[irLen,numIRs] = size(hL);
fVec = getFreqVec(Fs,irLen);
wVec = 2*pi*fVec;
nyqIndx = ceil((irLen+1)/2);

if find(strcmpi(varargin,'linearfit'))
    fitFlag = true;
    if nargin < 5
        fL = fVec(1);
        fU = fVec(nyqIndx);
    else
        if isempty(varargin{1})
            fL = fVec(1);
            fU = fVec(nyqIndx);
        else
            validateattributes(varargin{1},{'numeric'},{'vector',...
                'nonnan','finite','real','nonnegative','<=',...
                fVec(nyqIndx),'numel',2},'estimateIPD','[FL,FU]',4);
            fL = varargin{1}(1);
            fU = varargin{1}(2);
        end
    end
else
    fitFlag = false;
    if nargin < 4
        fL = fVec(1);
        fU = fVec(nyqIndx);
    else
        validateattributes(varargin{1},{'numeric'},{'vector','nonnan',...
            'finite','real','nonnegative','<=',fVec(nyqIndx),'numel',2},...
            'estimateIPD','[FL,FU]',4);
        fL = varargin{1}(1);
        fU = varargin{1}(2);
    end
end

% === COMPUTE IPD === %

pL = getPhaseSpec(hL,'unwrap');
pR = getPhaseSpec(hR,'unwrap');
dP = pR-pL;
IPD = zeros(nyqIndx,numIRs);
IPD(2:nyqIndx,:) = dP(2:nyqIndx,:)./repmat(wVec(2:nyqIndx),1,numIRs);
if tpFlag
    varargout{1} = IPD.';
else
    varargout{1} = IPD;
end

% === COMPUTE ITD IF REQUESTED === %

if nargout > 1
    if fitFlag
        ITD = estimateITD(hL,hR,Fs,'phase','range',[fL,fU],'linearfit');
    else
        ITD = estimateITD(hL,hR,Fs,'phase','range',[fL,fU]);
    end
    varargout{2} = ITD;
end

end
