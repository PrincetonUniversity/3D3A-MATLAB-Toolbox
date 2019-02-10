function onset = approximateHRIROnset(irIn,fS,varargin)
%APPROXIMATEHRIRONSET Approximate onset of head-related impulse responses.
%   B = APPROXIMATEHRIRONSET(A,FS) computes approximate onsets of the head-
%   related impulse responses (HRIRs) in A. If A is matrix, the individual 
%   HRIRs must be stored as column vectors. The sampling rate, FS, must 
%   also be specified in Hz. B is either a scalar or a row vector 
%   containing onsets specified in samples. The algorithm used is as
%   follows:
%       1. Window HRIRs using an asymmetric raised-cosine window with
%       approximately 1 ms pre-peak delay and post-peak duration.
%       2. Low-pass filter windowed HRIRs using minimum-phase, fourth order
%       butterworth filter with 2 kHz cut-off frequency.
%       3. Estimate onset as the sample corresponding to the maximum of the
%       cross-correlation between the filtered IR and its minimum phase
%       version.
%
%   B = APPROXIMATEHRIRONSET(...,'winParams',WP) specifies optional window
%   parameters, WP, as a 1-by-N cell array. Each element of the cell array
%   corresponds to a different window parameter. The first parameter must
%   be the window length in samples. The remaining parameters correspond to 
%   the optional Name-Value pair inputs of the windowSignal function and 
%   may be specified in any order.
%
%   B = APPROXIMATEHRIRONSET(...,'filterParams',FP) specifies optional
%   filter parameters, FP, as a 1-by-M cell array. The following may be
%   specified as Name-Value pairs within the cell array:
%       1. 'n',N - butterworth filter order
%       2. 'cutoff',fC - cut-off frequency in Hz.
%       3. 'fType',FT - type of filter ('low' - low-pass; 'high'-
%       high-pass; 'bandpass' - band-pass; 'stop' - band-stop)
%
%   B = APPROXIMATEHRIRONSET(...,'quiet','on') optionally specifies a flag 
%   that suppresses any warnings. By default, warnings are enabled.

%   EXAMPLE: For default behavior, the function call: 
%       B = APPROXIMATEHRIRONSET(A,FS)
%   is equivalent to:
%       B = APPROXIMATEHRIRONSET(A,FS,'winParams',{round(2*FS/1000),...
%           'wType',{'rc',[0,0.1]},'start',wS,'offset',oVec},...
%           'filterParams',{'n',4,'fType','low','cutoff',2000});
%   where wS and oVec are pre-computed.
%
%   See also ESTIMATEIRONSET.

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

narginchk(2,7);

% Parse and verify inputs
[inputs,extra] = parseApproximateHRIROnsetInputs(irIn,fS,varargin);

% Extract parsed inputs
irIn = inputs.irIn;
fS = inputs.fS;
WP = inputs.winParams;
FP = inputs.filterParams;
quietFlag = inputs.quiet;

% Check for unusually low value of fS
if fS < 1000 && strcmpi(quietFlag,'off')
    warning(['Detected unusually low input sampling rate. Note that',...
        ' sampling rate must be specified in Hz.'])
end

% Verify individual window parameters
validateattributes(WP{1,1},{'double'},{'scalar','nonempty','nonnan',...
    'finite','integer','positive'},'computeHRIROnset',...
    'first element in WP');

indx = find(strcmpi(WP,'wType'),1);
if indx && numel(WP) > indx
    validateattributes(WP{1,indx+1},{'cell'},{'nonempty','size',...
        [1,NaN]},'computeHRIROnset','input parameter for ''wType'' in WP');
    WPwType = WP{1,indx+1};
else
    WPwType = {'rc',[0,0.1]};
end

indx = find(strcmpi(WP,'start'),1);
if indx && numel(WP) > indx
    validateattributes(WP{1,indx+1},{'double'},{'scalar','nonempty',...
        'nonnan','finite','integer','positive'},'computeHRIROnset',...
        'input parameter for ''start'' in WP');
    WPstart = WP{1,indx+1};
else
    WPstart = extra{1,1};
end

indx = find(strcmpi(WP,'oVec'),1);
if indx && numel(WP) > indx
    validateattributes(WP{1,indx+1},{'double'},{'vector','nonempty',...
        'nonnan','finite','integer','nonnegative'},'computeHRIROnset',...
        'input parameter for ''oVec'' in WP');
    WPoVec = WP{1,indx+1};
else
    WPoVec = extra{1,2};
end

% Verify individual filter parameters
indx = find(strcmpi(FP,'n'),1);
if indx && numel(FP) > indx
    validateattributes(FP{1,indx+1},{'double'},{'scalar','nonempty',...
        'nonnan','finite','integer','positive'},'computeHRIROnset',...
        'input parameter for ''n'' in FP');
    FPn = FP{1,indx+1};
else
    FPn = extra{1,3};
end

indx = find(strcmpi(FP,'fType'),1);
if indx && numel(FP) > indx
    validateattributes(FP{1,indx+1},{'char'},{'scalartext','nonempty'},...
        'computeHRIROnset','input parameter for ''fType'' in FP');
    FPfType = FP{1,indx+1};
else
    FPfType = extra{1,5};
end

indx = find(strcmpi(FP,'cutoff'),1);
if indx && numel(FP) > indx
    validateattributes(FP{1,indx+1},{'double'},{'vector','nonempty',...
        'nonnan','finite','real','positive'},'computeHRIROnset',...
        'input parameter for ''cutoff'' in FP');
    FPcutoff = FP{1,indx+1};
    
    if strcmpi(FPfType,'bandpass') || strcmpi(FPfType,'stop')
        if length(FPcutoff) < 2 
            error(['When ''fType'' in FP is specified as ''bandpass''',...
                ' or ''stop'', the value for ''cutoff'' must be ',...
                'specified as [fL,fU] where fL < fU.'])
        else
            if any(FPcutoff(:) < 100) && strcmpi(quietFlag,'off')
                warning(['Detected low input filter cutoff frequency. ',...
                    'Note that cutoff frequency must be specified in Hz.'])
            end
            
            if length(FPcutoff) > 2
                error('The value for ''cutoff'' in FP must have length 2.')
            else
                if FPcutoff(1) >= FPcutoff(2)
                    error(['The value for ''cutoff'' in FP must be ',...
                        'specified as [fL,fU] where fL < fU.']);
                end
            end
        end
    else
        if any(FPcutoff(:) < 100) && strcmpi(quietFlag,'off')
            warning(['Detected low input filter cutoff frequency. Note',...
                ' that cutoff frequency must be specified in Hz.'])
        end
    end
else
    FPcutoff = extra{1,4};
end

% Main calculation

% 0. Compute irIn dimensions and generate impulse signal.
[irLen,~] = size(irIn);
imp = [1;zeros(irLen-1,1)];

% 1. Window HRIRs
[winIR,~,~] = windowSignal(irIn,WP{1,1},'wType',WPwType,'start',WPstart,...
    'offset',WPoVec);

% 2. Compute and apply minimum-phase low-pass filter
[num,den] = butter(FPn,FPcutoff/(fS/2),FPfType);
lpf_fir = filter(num,den,imp);
lpf_mpfir = makeMinPhaseIR(lpf_fir,'hilb');
winIR_lpf = filter(lpf_mpfir,1,winIR);

% 3. Estimate onset
onset = estimateIROnset(winIR_lpf,{'mpxc'});

end

function [inputs,extra] = parseApproximateHRIROnsetInputs(irIn,fS,opts)
%PARSEAPPROXIMATEHRIRONSETINPUTS Parse and verify inputs for the 
%approximateHRIROnset function.

p = inputParser;

% Required inputs
addRequired(p,'irIn',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite','real'},'approximateHRIROnset','A',1));
addRequired(p,'fS',@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','nonnan','finite','real','positive'},...
    'approximateHRIROnset','FS',2));

% If irIn is a vector, force it to be a column vector.
irIn = shiftdim(irIn);

% Compute default values for 'winParams'
wLen = round(2*fS/1000); % Approx. 2 ms.
maxOnsets = estimateIROnset(irIn,{'threshold',100});
minOnset = min(maxOnsets);
prePeakDelay = round(1*fS/1000); % Approx. 1 ms.
wS = max([1,minOnset-prePeakDelay+1]);
oVec = maxOnsets-minOnset;

% Specify default values for 'filterParams'
n = 4;
fType = 'low';
cutoff = 2000;

% Optional inputs
addParameter(p,'winParams',{wLen,'wType',{'rc',[0,0.1]},'start',wS,...
    'offset',oVec},@(x)validateattributes(x,{'cell'},{'nonempty','size',...
    [1,NaN]},'approximateHRIROnset','winParams'));
addParameter(p,'filterParams',{'n',n,'cutoff',cutoff,'fType',fType},...
    @(x)validateattributes(x,{'cell'},{'nonempty','size',[1,NaN]},...
    'approximateHRIROnset','filterParams'));
addParameter(p,'quiet','off',@(x)validateattributes(x,{'char'},{...
    'scalartext','nonempty'},'approximateHRIROnset',...
    'value for ''quiet'''));

% Specify extra variables to return
extra = {wS,oVec,n,cutoff,fType};

p.CaseSensitive = false;
p.FunctionName = 'approximateHRIROnset';

parse(p,irIn,fS,opts{:});

inputs = p.Results;

end
