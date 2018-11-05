function [sd,mSD] = computeSD(ref,test,DIM)
%COMPUTESD Compute spectral distortion
%   [sd,mSD] = COMPUTESD(ref,test) computes the spectral distortion between 
%   ref and test. Some possible options for ref and test might be:
%       1. frequency responses
%       2. magnitude responses
%   Both ref and test may be vectors, or matrices with the same dimensions.
%   The mean spectra distortion (mSD) is also returned. Both sd and mSD are 
%   returned in dB.
%
%   ___ = COMPUTESD(...,DIM) specifies the dimensions along which to
%   perform computation of mSD. See the SUM function for more information.

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
%       [1]. Xie (2013) - Head-Related Transfer Function and Virtual
%       Auditory Display, Second Edition, J. Ross. 

narginchk(2,3);

if nargin < 3
    DIM = 1;
end

validateattributes(ref,{'double'},{'2d','nonempty','nonnan','finite',...
    'ndims',2},'computeSD','ref',1);
validateattributes(test,{'double'},{'2d','nonempty','nonnan','finite',...
    'size',size(ref)},'computeSD','test',2);
validateattributes(DIM,{'double'},{'scalar','nonempty','nonnan',...
    'finite','integer','positive','<=',2},'computeSD','DIM',3);

ref = shiftdim(ref); % If row vector, convert to column vector.
test = shiftdim(test);

% The following formulas were adapted from page 138 of Xie [1].
sd = mag2db(abs(ref)./abs(test));
mSD = rms(sd,DIM);

end
