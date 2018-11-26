function [az,ITD,ILD,fVec] = wrap_may2011(bIn,fS,varargin)
%WRAP_MAY2011 Wrapper for may2011 from the Auditory Modeling Toolbox.
%   [az,ITD,ILD,fVec] = WRAP_MAY2011(bIn,fS) takes a steady-state binaural 
%   signal, bIn, and estimates azimuth, ITD, and ILD in the frontal 
%   horizontal plane using the binaural localization model by May et al. 
%   [1], as implemented in the Auditory Modeling Toolbox (AMT). The 
%   sampling rate of the binaural signal must be specified in Hz. The 
%   output azimuth, az, is a structure consisting of:
%       1. the raw frequency-dependent azimuth estimate (see out.azimuth in
%       AMT v0.9.9 documentation of MAY2011.)
%       2. the mean of the above raw estimates
%   All azimuth values are specified in degrees and follow the conventions
%   of the CIPIC interaural coordinate system. The output ITD and ILD are 
%   also structures consisting of:
%       1. the raw frequency-dependent ITD/ILD estimates (see out.itd and
%       out.ild in AMT v0.9.9 documentation of MAY2011.)
%       2. the mean of the above raw estimates
%   All ITD values are specified in seconds, while ILD values are specified
%   in dB. Negative ITD and ILD corresponds to a sound source on the left 
%   (i.e. negative azimuth) - see page 3 of May et al. [1]. A vector of 
%   frequencies in Hz, fVec, at which the raw azimuth, ITD and ILD 
%   estimates are calculated is also provided. An alternative way to 
%   specify this command is: WRAP_MAY2011(bIn,fS,'bInType','ss').
%
%   ___ = WRAP_MAY2011(...,'bInType','imp') indicates that bIn consists of 
%   binaural impulse responses. In this case, steady-state binaural signals 
%   are generated internally by convolving each channel in bIn with a pink 
%   noise signal.
%
%   ___ = WRAP_MAY2011(...,'stimType',TYPE) indicates the type of stimulus
%   signal to use to generate steady-state binaural signals when 'bInType'
%   is specified as 'imp'. The following may be specified for TYPE:
%       1. 'pink-noise' - pink noise signal (default)
%       2. 'ess' - exponential sine sweep signal
%   This command is only applicable when 'bInType' is set to 'imp'.
%
%   Needs: Auditory Modeling Toolbox (AMT).
%
%   See also MAY2011.

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

%   References:
%       [1]. May et al. (2011) - A Probabilistic Model for Robust 
%       Localization Based on a Binaural Auditory Front-End.

narginchk(2,6);

if exist('amt_version','file') ~= 2
    error('The Auditory Modeling Toolbox must be in the MATLAB path.');
end

% The following check is added because the auto-download feature of the
% Auditory Modeling Toolbox (to download modeldata.mat automatically) does 
% not work reliably.
if exist('modeldata.mat','file') ~= 2
    error(['Requires modeldata.mat in /code/auxdata/may2011 within the',...
        ' Auditory Modeling Toolbox. Download from: ',...
        'http://www.sofacoustics.org/data/amt-0.9.9/auxdata/may2011/',...
        'modeldata.mat']);
end

% Parse and verify inputs
inputs = parseWrap_May2011Inputs(bIn,fS,varargin);

% Extract parsed inputs
bIn = inputs.bIn;
fS = inputs.fS;
bInType = inputs.bInType;
stimType = inputs.stimType;
inputLen = ceil(0.2*fS);
switch lower(bInType)
    case 'imp' % Generate steady-state binaural signal from input IRs
        switch lower(stimType)
            case 'pink-noise'
                pinkNoiseSource = dsp.ColoredNoise(1,inputLen,1);
                stimulus = step(pinkNoiseSource);
                stimulus = stimulus/max(abs(stimulus));
            case 'ess'
                stimulus = exponentialSineSweep(1,fS/2,fS,inputLen/fS,...
                    'type','phase');
            otherwise
                error('Invalid TYPE specification for option: stimType.')
        end
        
        input = fftConv(bIn,stimulus,'lin');
    case 'ss' % Assume input binaural signal is steady-state        
        input = [bIn; zeros(inputLen,2)];
    otherwise
        error(['Unrecognized third input. Only ''ss'' and',...
            ' ''imp'' are accepted.']);
end

maxInput = max(max(abs(input)));
input = input/maxInput;
out = may2011(input,fS); % Call to function in AMT.
fVec = erbspace(80,5000,32).'; % See page 2 of May et al. [1].
az.raw = mean(out.azimuth,2,'omitnan'); % Frequency-dependent azimuth.
az.mean = logmean(az.raw,fVec,[fVec(2),fVec(10)]);
ITD.raw = out.itd; % Frequency-dependent ITD estimates.
ITD.mean = mean(out.itd);
ILD.raw = out.ild; % Frequency-dependent ILD estimates.
ILD.mean = mean(out.ild);

end

function inputs = parseWrap_May2011Inputs(bIn,fS,opts)
%PARSEWRAP_MAY2011INPUTS Parse and verify inputs for the wrap_may2011 
%function.

p = inputParser;

% Required inputs
addRequired(p,'bIn',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite','real','size',[NaN,2]},'wrap_may2011',...
    'bIn',1));
addRequired(p,'fS',@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','positive','real'},'wrap_may2011','fS',2));

% Optional inputs
addParameter(p,'bInType','ss',@(x)validateattributes(x,{'char'},...
    {'scalartext'},'wrap_may2011','type of input signal'));
addParameter(p,'stimType','pink-noise',@(x)validateattributes(x,...
    {'char'},{'scalartext'},'wrap_may2011','type of stimulus signal'));

p.CaseSensitive = false;
p.FunctionName = 'wrap_may2011';

parse(p,bIn,fS,opts{:});

inputs = p.Results;

end
