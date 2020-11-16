function ITD = estimateITD(hL,hR,Fs,METHOD,varargin)
%ESTIMATEITD Estimate the interaural time difference (ITD), in seconds.
%   ITD = ESTIMATEITD(hL,hR,Fs) estimates ITD for input HRIRs hL and hR at 
%   sampling rate Fs specified in Hz using thresholding with a threshold of
%   0.2 (20% of the absolute maximum of each IR). Negative ITDs correspond
%   to sound sources on the left. hL and hR may be vectors or matrices but
%   must have the same dimensions. If matrices, the IRs must be stored as
%   columns.
%
%   ITD = ESTIMATEITD(...,METHOD) additionally specifies the method of 
%   estimating ITD. The options are:
%       1. 'group delay' or 'group' - estimate ITD accurate to 1 sample at 
%   the specified sampling rate by computing the difference in the group
%   delays of hL and hR, and averaging the result over nominal frequencies 
%   of 0 to 1500 Hz. To modify the range over which averaging is performed,
%   see the additional options below.
%       2. 'cross-correlation' or 'xcorr' estimates ITD accurate to 1 
%   sample at the specified sampling rate by determining the first sample 
%   position at which the cross-correlation between hL and hR reaches an 
%   absolute maximum.
%       3. 'phase delay' or 'phase' estimates ITD with fractional-sample 
%   accuracy by determing the difference between the unwrapped phase at the
%   right ear and that at the left ear, each averaged over 0 to 1500 Hz
%   nominal frequencies. To modify the range over which averaging is 
%   performed, see the additional options below.
%       4. 'thresholding' or 'thresh' estimates ITD accurate to 1 sample at 
%   the specified sampling rate by thresholding with a threshold of 0.2 
%   (20% of the absolute maximum of each IR). This is the default.
%       5. 'mpxc' estimates ITD accurate to 1 sample at the specified
%   sampling rate by first computing the onset of hL and hR by 
%   cross-correlating each with its corresponding minimum-phase version, 
%   and then taking the difference between these computed onsets. For more 
%   on the cross-correlation algorithm, see ESTIMATEIRONSET.
%       6. 'mpxc_robust' estimates ITD accurate to a fraction of a sample
%   by first computing the onset of windowed versions of hL and hR by 
%   cross-correlating each with its corresponding minimum-phase version, 
%   and then taking the difference between these computed onsets. For more 
%   on the cross-correlation algorithm, see ESTIMATEIRONSET.
%
%   ITD = ESTIMATEITD(...,Name1,Value1,...) specifies optional 
%   comma-separated pairs of Name,Value arguments, where Name is the 
%   argument name and Value is the corresponding value. Name must appear 
%   inside single quotes (' '). You can specify the following:
%   
%   'filter'        Specifications for a Butterworth filter to filter hL
%                   and hR prior to ITD estimation. The specifications must
%                   be a 1x3 cell array specifying {filter order, cut-off
%                   frequency in Hz, filter type} for a Butterworth
%                   filter. For order and type specifications, see the
%                   documentation for 'butter'. The three parameters must
%                   be specified in the order shown above. An optional
%                   fourth parameter may be specified as either 'zerophase'
%                   to design a zero-phase filter, or 'minphase' to design
%                   a minimum-phase filter.
%
%   'resample'      Resample hL and hR prior to computing ITD. This must
%                   be a scalar > 0 such that a value in the range (0,1)
%                   corresponds to downsampling, and a value > 1
%                   corresponds to upsampling.
%
%   'range'         Frequency range, in Hz, over which averaging is
%                   performed for the 'group delay' and 'phase delay' 
%                   estimation methods. This must be specified as a 1x2 or
%                   a 2x1 vector with the first value specifying the lower
%                   bound.
%
%   'threshold'     Threshold fraction used for the thresholding approach
%                   to estimate ITD. This must be a scalar between 0 and 1,
%                   with 1 corresponding to the maximum absolute value of a
%                   given IR, and 0 corresponding to the first non-zero
%                   sample.
%
%   'linearfit'     [No Value] Specifying this parameter first generates a
%                   linear fit to the unwrapped phase response over the 
%                   frequency range specified in 'range' (or the default of
%                   0 to 1500 Hz) prior to estimating ITD using the 'phase 
%                   delay' approach.
%
%   Needs: Signal Processing Toolbox.
%
%   See also ESTIMATEIRONSET, THRESHOLDIRS.

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

narginchk(3,13);

if nargin < 4
    METHOD = 'thresholding';
end

validateattributes(hL,{'numeric'},{'2d','nonempty','nonnan','finite'},...
    'estimateITD','hL',1);
validateattributes(hR,{'numeric'},{'2d','nonempty','nonnan','finite',...
    'size',size(hL)},'estimateITD','hR',2);
validateattributes(Fs,{'numeric'},{'scalar','nonempty','nonnan',...
    'finite','real','positive'},'estimateITD','Fs',3);

hL = shiftdim(hL);
hR = shiftdim(hR);

% Perform optional filtering first
indx = find(strcmpi(varargin,'filter'),1);
if ~isempty(indx)
    specCell = varargin{indx+1};
    n = specCell{1,1}; % filter order; see 'butter' help
    cutoff = specCell{1,2}; % cutoff frequency in Hz
    ftype = specCell{1,3}; % filter type; see 'butter' help
    Wn = cutoff/(Fs/2);
    [z,p,k] = butter(n,Wn,ftype); % From Signal Processing Toolbox
    [b,a] = zp2tf(z,p,k);
    if length(specCell) == 4 && strcmpi(specCell{1,4},'zerophase')
        hL = filtfilt(b,a,hL); % From Signal Processing Toolbox
        hR = filtfilt(b,a,hR); % From Signal Processing Toolbox
    elseif length(specCell) == 4 && strcmpi(specCell{1,4},'minphase')
        irLen = size(hL,1);
        imp = [1;zeros(irLen-1,1)];
        lpf_imp = filter(b,a,imp);
        lpf_imp_mp = makeMinPhaseIR(lpf_imp,'hilb');
        hL = filter(lpf_imp_mp,1,hL);
        hR = filter(lpf_imp_mp,1,hR);
    else
        hL = filter(b,a,hL);
        hR = filter(b,a,hR);
    end
    filtFlag = true;
else
    filtFlag = false;
end

% Perform resampling next ('upsample' included for backwards-compatibility)
indx = find(strcmpi(varargin,'resample') | strcmpi(varargin,'upsample'),1);
if ~isempty(indx)
    [p,q] = rat(varargin{indx+1});
    [~,aaf] = resample([hL,hR],p,q); % From Signal Processing Toolbox
    if strcmpi(METHOD,'mpxc')
        aaf_mp = makeMinPhaseIR(aaf,'hilb');
        hL = resample(hL,p,q,aaf_mp);
        hR = resample(hR,p,q,aaf_mp);
    else
        hL = resample(hL,p,q,aaf);
        hR = resample(hR,p,q,aaf);
    end
    Fs = Fs*(p/q);
end

indx = find(strcmpi(varargin,'range'),1);
if ~isempty(indx)
    if numel(varargin{indx+1})==1
        fU = varargin{indx+1};
        fL = 0;
    else
        fL = varargin{indx+1}(1);
        fU = varargin{indx+1}(2);
    end
else
    fL = 0;
    fU = 1500;
end

indx = find(strcmpi(varargin,'threshold'),1);
if ~isempty(indx)
    thp = varargin{indx+1};
else
    thp = 0.2;
end

switch lower(METHOD)
    case {'group delay','group'}
        dL = estimateIROnset(hL,{'grpdelay',Fs,[fL,fU]});
        dR = estimateIROnset(hR,{'grpdelay',Fs,[fL,fU]});
        d = dL-dR;
    case {'cross-correlation','xcorr'}
        numIRs = size(hL,2);
        d = zeros(1,numIRs);
        for ii = 1:numIRs
            [xc,lagVec] = xcorr(hL(:,ii),hR(:,ii));
            [~,lagIndex] = max(abs(xc));
            d(ii) = lagVec(lagIndex);
        end
    case {'phase delay','phase'}
        indx = find(strcmpi(varargin,'linearfit'),1);
        if ~isempty(indx)
            dL = estimateIROnset(hL,{'phase',Fs,[fL,fU],'linearfit'});
            dR = estimateIROnset(hR,{'phase',Fs,[fL,fU],'linearfit'});
        else
            dL = estimateIROnset(hL,{'phase',Fs,[fL,fU]});
            dR = estimateIROnset(hR,{'phase',Fs,[fL,fU]});
        end
        d = dR-dL;
    case {'thresholding','thresh'}
        dL = estimateIROnset(hL,{'threshold',thp*100});
        dR = estimateIROnset(hR,{'threshold',thp*100});
        d = dL-dR;
    case 'mpxc'
        dL = estimateIROnset(hL,{'mpxc'});
        dR = estimateIROnset(hR,{'mpxc'});
        d = dL-dR;
    case 'mpxc_robust'
        if ~filtFlag
            d = estimateITD(hL,hR,Fs,'mpxc_robust','filter',{4,2000,...
                'low'});
            d = d*Fs;
        else
            dL = estimateIROnset(hL,{'mpxc_robust'});
            dR = estimateIROnset(hR,{'mpxc_robust'});
            d = dL-dR;
        end
    otherwise
        error('Invalid input for METHOD.')
end

ITD = d/Fs;

end
