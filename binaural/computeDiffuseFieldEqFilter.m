function fIR = computeDiffuseFieldEqFilter(h,fS,varargin)
%COMPUTEDIFFUSEFIELDEQFILTER Compute diffuse-field equalization filter.
%   fIR = COMPUTEDIFFUSEFIELDEQFILTER(h,fS) computes a diffuse-field
%   equalization filter for a set of diffuse-field IRs given in h. h must 
%   be a matrix with columns containing IRs corresponding to different
%   directions. The sampling rate, fS, must be specified in Hz. The output,
%   fIR, is a vector of length equal to that of the IRs in h. To compute
%   diffuse-field IRs required as input to this function, see 
%   COMPUTEDIFFUSEFIELDIR.
%
%   fIR = COMPUTEDIFFUSEFIELDEQFILTER(...,FILTERLEN) optionally specifies
%   the desired length of fIR in samples. Specify FILTERLEN as [] to use
%   the default value of size(h,1) when also specifying other options (see
%   below).
%
%   fIR = COMPUTEDIFFUSEFIELDEQFILTER(...,'avgFreqRange',[FL,FU]) specifies
%   the frequency range over which avg. spectrum level (required for
%   computing the inverse filter) is estimated. Default: [200,18000]
%
%   fIR = COMPUTEDIFFUSEFIELDEQFILTER(...,'invFreqRange',[FL,FU]) specifies
%   the frequency range over which inversion is performed.
%   Default: [20,20000]
%
%   See also COMPUTEDIFFUSEFIELDIR.

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
%   Copyright (c) 2021 Princeton University
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

if nargin < 3
    FILTERLEN = size(h,1); % Default fIR length in samples.
else
    FILTERLEN = varargin{1};
    if isempty(FILTERLEN)
        FILTERLEN = size(h,1); % Default fIR length in samples.
    end
end

hLen = size(h,1);
invFiltLen = 2*max([hLen,FILTERLEN,ceil(fS/20)])+1;
if invFiltLen > hLen
    hPad = padarray(h,[invFiltLen-hLen,0],0,'post');
end

indx = find(strcmpi(varargin,'avgFreqRange'),1);
if isempty(indx)
    w1a = 200/(fS/2);
    w2a = 18000/(fS/2);
else
    w1a = min(varargin{indx+1})/(fS/2);
    w2a = max(varargin{indx+1})/(fS/2);
end

indx = find(strcmpi(varargin,'invFreqRange'),1);
if isempty(indx)
    w1i = 20/(fS/2);
    w2i = 20000/(fS/2);
else
    w1i = min(varargin{indx+1})/(fS/2);
    w2i = max(varargin{indx+1})/(fS/2);
end

eqFilterIR = computeInverseFilter(hPad,'custom',{'avgFreqRange',...
    [w1a,w2a],'invFreqRange',[w1i,w2i],'maxDynRange',24});
eqFilterIR = shiftSignal(ifft(abs(fft(eqFilterIR)),'symmetric'),...
    round(FILTERLEN/2));
fIR = windowSignal(eqFilterIR,FILTERLEN,'wType',{'tukey',0.1});

end
