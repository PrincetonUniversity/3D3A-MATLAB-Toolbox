function hMat_diffuseField = computeDiffuseFieldIR(hMat,METHOD)
%COMPUTEDIFFUSEFIELDIR Compute diffuse-field impulse responses.
%   hMat_diffuseField = COMPUTEDIFFUSEFIELDIR(hMat) computes the
%   diffuse-field impulse responses (IRs) from the input IRs in hMat by
%   averaging the individual power spectra. hMat must contain IRs stored as
%   columns. The computed diffuse-field IRs have a linear phase response
%   with the maximum of the absolute magnitude occurring at the sample
%   value corresponding to half the length of the IRs.
%
%   hMat_diffuseField = COMPUTEDIFFUSEFIELDIR(hMat,METHOD) optionally
%   specifies the method by which to compute the diffuse-field IRs. METHOD
%   must be a cell array with the first element being the MethodName and
%   following elements being parameters specific to the chosen method. The
%   following may be specified for METHOD:
%       1. {'avgPS'} - unweighted average of individual power spectra.
%       2. {'avgMagdB'} - unweighted average of individual magnitude
%       spectra in dB.
%
%   See also COMPUTEDIFFUSEFIELDEQFILTER.

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

narginchk(1,2);

if nargin < 2
    METHOD = {'avgPS'};
end

hMat = shiftdim(hMat);
switch lower(METHOD{1})
    case 'avgps'
        psMat = (getMagSpec(hMat)).^2;
        magMat_diffuseField = sqrt(mean(psMat,2,'omitnan'));
    case 'avgmagdb'
        magdBMat = getMagSpecdB(hMat);
        magMat_diffuseField = db2mag(mean(magdBMat,2,'omitnan'));
    otherwise
        error('Invalid METHOD specification.');
end

irLen = size(magMat_diffuseField,1);
hMat_diffuseField = ifft(magMat_diffuseField,'symmetric');
% Make output causal
hMat_diffuseField = shiftSignal(hMat_diffuseField,irLen/2);

end
