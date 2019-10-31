function [test_lm,gain_dB] = levelMatchSignal(ref,test,Fs,STIM)
%LEVELMATCHSIGNAL Match the level of two signals.
%   [P,G] = LEVELMATCHSIGNAL(R,T,FS) matches the avg. level of the signal
%   in T to that of R and returns the re-normalized version of T as P. Also
%   returned is the dynamic range loss, G, in dB, due to the re-
%   normalization. Average levels are computed by convolving R/T with pink 
%   noise (to avoid this step, see options below), filtering the resulting 
%   signal through critical band filters with center frequencies spaced 1 
%   ERB apart, and finally computing the average level in those critical 
%   bands whose center frequencies are between 1 and 5 kHz. R and T must
%   both be vectors. P is a vector and G is a scalar. The sampling rate,
%   FS, must be specified in Hz.
%
%   [P,G] = LEVELMATCHSIGNAL(R,T,FS,STIM) optionally specifies STIM which
%   can take the following options:
%       (i) 'pink' (default) - R and T are convolved with pink noise prior 
%       to level matching. Typically, this means R and T are impulse
%       responses.
%       (ii) 'none' - R and T are used as is (no convolution is performed).
%       This option is suitable if R and T are to be treated as signals
%       instead of impulse responses.

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
%   Copyright (c) 2019 Princeton University
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

narginchk(3,4);

if nargin < 4
    STIM = 'pink';
end

% Validate format of inputs
validateattributes(ref,{'double'},{'vector','nonempty','nonnan',...
    'finite','real'},'levelMatchSignal','R',1)
validateattributes(test,{'double'},{'vector','nonempty','nonnan',...
    'finite','real'},'levelMatchSignal','T',2)
validateattributes(Fs,{'double'},{'scalar','nonempty','nonnan',...
    'finite','real','positive'},'levelMatchSignal','Fs',3)
validateattributes(STIM,{'char'},{'scalartext','nonempty'},...
    'levelMatchSignal','STIM',4)
ref = shiftdim(ref);
test = shiftdim(test);

switch lower(STIM)
    case 'pink'
        % Set stimulus characteristics and generate raw stimulus
        
        T = 1; % Duration in seconds.
        rawStim = generateStimulus('pink',Fs,T); % Raw pink noise
        rawStimLen = length(rawStim);
        winStim = windowSignal(rawStim,rawStimLen,'wType',{'tukey',0.2});
        refSig = fftConv(ref,winStim,'lin');
        testSig = fftConv(test,winStim,'lin');
    case 'none'
        refSig = ref;
        testSig = test;
    otherwise
        error('Unrecognized input for STIM.')
end

% Generate critical band filters
erbFCs = getERBFreqVec(1000,5000);
impVec = zeros(2^nextpow2(0.04*Fs),1);
impVec(1) = 1;
cbFiltMat = applyCBFilter(impVec,Fs,erbFCs);
filtLen = size(cbFiltMat,1);
cbFiltMat = shiftSignal(cbFiltMat,filtLen/2);

% Filter ref and test signals using critical band filters
filtRef = fftConv(cbFiltMat,refSig,'lin');
filtTest = fftConv(cbFiltMat,testSig,'lin');
filtRefLen = size(filtRef,1);
filtTestLen = size(filtTest,1);

% Compute average perceived power 
energyRef = mean(sum(abs(filtRef).^2))/filtRefLen;
energyTest = mean(sum(abs(filtTest).^2))/filtTestLen;

% Computing level matching gain and apply to test signal
gain = sqrt(energyRef/energyTest);
gain_dB = mag2db(gain);
test_lm = test*gain;

end
