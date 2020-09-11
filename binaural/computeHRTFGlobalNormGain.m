function G = computeHRTFGlobalNormGain(HL,HR,S,Fs,varargin)
%COMPUTEHRTFGLOBALNORMGAIN Estimate global normalization gain for HRTF
%dataset.
%   G = COMPUTEHRTFGLOBALNORMGAIN(HL,HR,S,Fs) takes a full HRTF dataset
%   consisting of left-ear HRTFs, HL, and right-ear HRTFs, HR along with a
%   source position matrix, S, and computes a global normalization gain, G,
%   to apply to the HRTFs in HL and HR so that the estimated average 
%   magnitude (up to 1 kHz) is 0 dB for a response that might have been 
%   measured at the center of the head (with the head absent). Such a
%   response is estimated by averaging the left- and right-ear HRTFs for 
%   the frontal direction (i.e., Az = 0 and El = 0 in SOFA spherical 
%   coordinates), or the closest available such direction.
%
%       HL and HR must each be an N-by-P matrix where N is the length of an
%       individual HRIR and P is the number of source directions. The data
%       in H must consist of impulse responses.
%
%       S must be a P-by-3 matrix of source positions specified in SOFA
%       cartesian coordinates.
%
%       Fs is the sampling rate in Hz.
%
%   The output, G, will be a scalar specifying the gain, in dB, to apply to 
%   both HL and HR.
%
%   G = COMPUTEHRTFGLOBALNORMGAIN(___,'fc',Fc) specifies the frequency,
%   Fc, in Hz such that the average magnitude estimate described above is 
%   computed from 0 to (approximately) Fc Hz. As described above, the 
%   default value of Fc is 1 kHz.
%
%   Note on the algorithm: This approach to normalization is used because
%   the magnitude response of the rigid-sphere HRTF with antipodal ears is
%   very close to 0 dB for frontal sources, *irrespective of source
%   distance*. When the source is infinitely far away, the low-frequency
%   magnitude response is close to 0 dB irrespective of source direction. 
%   But for non-frontal sources and finite source distances (typically the 
%   case for acoustically-measured HRTFs), low-frequency magnitudes can
%   depart from 0 dB quite significantly (depending on source distance).

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

narginchk(4,6);

% Validate inputs
validateattributes(HL,{'numeric'},{'2d','nonempty','finite','nonnan'},...
    'computeHRTFGlobalNormGain','HL',1);
validateattributes(HR,{'numeric'},{'2d','nonempty','finite','nonnan'},...
    'computeHRTFGlobalNormGain','HR',2);
HL = shiftdim(HL); % If HL is a vector, force it to be a column vector.
HR = shiftdim(HR);
if size(HL) ~= size(HR)
    error(['HL and HR have different sizes. HL and HR must have the',...
        ' same size.'])
end
[HLen,numPos] = size(HL);
validateattributes(S,{'numeric'},{'2d','nonempty','finite','nonnan',...
    'real','size',[numPos,3]},'computeHRTFGlobalNormGain','S',3);
validateattributes(Fs,{'numeric'},{'scalar','nonempty','finite','real',...
    'nonnan','positive'},'computeHRTFGlobalNormGain','Fs',4);

% Parse optional inputs
fcIndx = find(strcmpi(varargin,'fc'),1);
if ~isempty(fcIndx)
    Fc = varargin{fcIndx+1};
    validateattributes(Fc,{'double'},{'scalar','nonempty','finite',...
        'nonnan','real','nonnegative'},'computeHRTFGlobalNormGain','Fc');
else
    Fc = 1000; % In Hz
end

% Extract HRIRs corresponding to frontal position
[HL_front,~,front_indx] = extractIR(HL,S,[1,0,0],true);
HR_front = HR(:,front_indx);

% Compute G
magRespL = getMagSpec(HL_front);
magRespR = getMagSpec(HR_front);
magRespMean = mean([magRespL,magRespR],2);
fVec = getFreqVec(Fs,HLen);
[~,f_cutoff_indx] = min(abs(fVec-Fc));
G = -mag2db(mean(magRespMean(1:f_cutoff_indx)));

end
