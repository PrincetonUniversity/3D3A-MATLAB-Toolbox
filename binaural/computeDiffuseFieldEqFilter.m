function fIR = computeDiffuseFieldEqFilter(h,fS,FILTERLEN)
%COMPUTEDIFFUSEFIELDEQFILTER Compute diffuse-field equalization filter.
%   fIR = COMPUTEDIFFUSEFIELDEQFILTER(h,fS) computes a diffuse-field
%   equalization filter for a set of diffuse-field IRs given in h. h must 
%   be a matrix with columns containing IRs corresponding to different
%   directions. The sampling rate, fS, must be specified in Hz. The output,
%   fIR, is a vector of length 5 ms (at the specified fS). To compute
%   diffuse-field IRs required as input to this function, see 
%   COMPUTEDIFFUSEFIELDIR.
%
%   fIR = COMPUTEDIFFUSEFIELDEQFILTER(...,FILTERLEN) optionally specifies
%   the desired length of fIR in samples.
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

narginchk(2,3);

if nargin < 3
    FILTERLEN = size(h,1); % Default fIR length in samples.
end

w1 = 500/(fS/2);
w2 = 15000/(fS/2);
eqFilterIR = computeInverseFilter(h,'gardner1994',{'avgRange',[w1,w2]});

fIR = makeMinPhaseIR(eqFilterIR,'hilb');
fIR = windowSignal(fIR,FILTERLEN,'wType',{'rc',[0,0.5]});

end
