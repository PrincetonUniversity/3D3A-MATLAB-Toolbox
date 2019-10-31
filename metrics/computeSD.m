function sD = computeSD(ref,test)
%COMPUTESD Compute spectral distortion
%   D = COMPUTESD(R,T) computes the distortion, in dB, between 
%   corresponding spectra in R and T. Both R and T may be vectors, or 
%   matrices with the same dimensions. R and T may contain either
%   complex-valued frequency responses or magnitude spectra. COMPUTESD then
%   computes, element-wise, the ratio |R|/|T| and returns the result in dB.

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

%   Ref:
%       [1]. Xie (2013) - Head-Related Transfer Function and Virtual
%       Auditory Display, Second Edition, J. Ross. 

narginchk(2,2);

% Check inputs
validateattributes(ref,{'double'},{'2d','nonempty','nonnan','finite',...
    'ndims',2},'computeSD','R',1);
validateattributes(test,{'double'},{'2d','nonempty','nonnan','finite',...
    'size',size(ref)},'computeSD','T',2);
ref = shiftdim(ref); % If row vector, convert to column vector.
test = shiftdim(test);

% The following formula was adapted from page 138 of Xie [1].
sD = mag2db(abs(ref)./abs(test));

end
