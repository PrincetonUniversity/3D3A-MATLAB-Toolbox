function [az,ITD,ILD,fVec] = wrap_may2011(bIn,fS,varargin)
%WRAP_MAY2011 Wrapper for may2011 from the Auditory Modeling Toolbox.
%   [az,ITD,ILD,fVec] = WRAP_MAY2011(bIn,fS) takes a binaural signal, bIn,
%   and estimates azimuth, ITD, and ILD in the frontal horizontal plane 
%   using the binaural localization model by May et al. [1], as implemented 
%   in the Auditory Modeling Toolbox (AMT). The sampling rate, fS, of the 
%   binaural signal must be specified in Hz. The output azimuth, az, is a 
%   matrix consisting of frequency-dependent azimuth estimates (see 
%   out.azimuth in the AMT documentation for MAY2011.), with each row 
%   corresponding to a different frequency, and each column to a different 
%   time-domain frame. All azimuth values are specified in degrees and 
%   follow the conventions of the CIPIC interaural coordinate system. 
%
%   The output ITD and ILD are also matrices (with the same structure as 
%   the az matrix) consisting of the raw frequency-dependent ITD/ILD 
%   estimates (see out.itd and out.ild in the AMT documentation for 
%   MAY2011.). All ITD values are specified in seconds, while ILD values 
%   are specified in dB. Negative ITD and ILD correspond to a sound source 
%   on the left (i.e. negative azimuth) - see page 3 of May et al. [1]. A 
%   vector of frequencies in Hz, fVec, at which the raw azimuth, ITD and 
%   ILD estimates are calculated is also provided. An alternative way to 
%   specify this command is: WRAP_MAY2011(bIn,fS,'bInType','sig').
%
%   ___ = WRAP_MAY2011(...,'bInType','imp') indicates that bIn consists of 
%   binaural impulse responses. In this case, binaural signals are 
%   generated internally by convolving each channel in bIn with a white 
%   noise signal.
%
%   ___ = WRAP_MAY2011(...,'stimType',TYPE) indicates the type of stimulus
%   signal to use to generate binaural signals when 'bInType' is specified 
%   as 'imp'. The following may be specified for TYPE:
%       1. 'white-noise' - white noise signal (default)
%       2. 'pink-noise' - pink noise signal
%       3. 'ess' - exponential sine sweep signal
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

if exist('amt_start','file') ~= 2
    error('The Auditory Modeling Toolbox must be in the MATLAB path.');
end

% Add auxdata folder in AMT-Toolbox to MATLAB path since amt_start does not
% do this.
mainPath = fileparts(which('amt_start.m'));
auxDataPath = fullfile(mainPath,'auxdata');
if exist(auxDataPath,'dir') == 7
    addpath(genpath(auxDataPath));
else
    mkdir(fullfile(auxDataPath,'may2011'));
end

% Parse and verify inputs
inputs = parseWrap_May2011Inputs(bIn,fS,varargin);

% Extract parsed inputs
bIn = inputs.bIn;
fS = inputs.fS;
bInType = inputs.bInType;
stimType = inputs.stimType;
inputLen = ceil(0.02*fS);
bInLen = size(bIn,1);
switch lower(bInType)
    case 'imp' % Generate steady-state binaural signal from input IRs
        switch lower(stimType)
            case 'white-noise'
                whiteNoiseSource = dsp.ColoredNoise(0,inputLen,1);
                stimulus = step(whiteNoiseSource);
                stimulus = stimulus/max(abs(stimulus));
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
    case 'sig' % Assume input binaural signal is not an impulse response       
        input = [bIn; zeros(inputLen-mod(bInLen,inputLen),2)];
    otherwise
        error(['Unrecognized third input. Only ''sig'' and',...
            ' ''imp'' are accepted.']);
end

maxInput = max(max(abs(input)));
input = input/maxInput;
out = may2011(input,fS); % Call to function in AMT.
fVec = erbspace(80,5000,32).'; % See page 2 of May et al. [1].
az = out.azimuth; % Frequency-dependent azimuth.
ITD = out.itd; % Frequency-dependent ITD estimates.
ILD = out.ild; % Frequency-dependent ILD estimates.

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
addParameter(p,'bInType','sig',@(x)validateattributes(x,{'char'},...
    {'scalartext'},'wrap_may2011','type of input signal'));
addParameter(p,'stimType','white-noise',@(x)validateattributes(x,...
    {'char'},{'scalartext'},'wrap_may2011','type of stimulus signal'));

p.CaseSensitive = false;
p.FunctionName = 'wrap_may2011';

parse(p,bIn,fS,opts{:});

inputs = p.Results;

end
