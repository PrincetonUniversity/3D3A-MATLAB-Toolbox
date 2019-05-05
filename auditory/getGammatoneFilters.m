function h = getGammatoneFilters(fc,Fs,IRLen,fType)
%GETGAMMATONEFILTERS Auditory filters that represent human critical bands.
%   H = GETGAMMATONEFILTERS(FC,FS,N) or H = GETGAMMATONEFILTERS(FC,FS,N,...
%   'complex') computes length N impulse responses of gammatone filters 
%   with center frequencies specified by the vector FC. By default, the 
%   impulse responses (IRs) are complex-valued. To generate real-valued 
%   IRs, see below. 
%
%   H = GETGAMMATONEFILTERS(FC,FS,N,'real') generates real-valued IRs
%   instead. The returned IRs are causal.
%   
%   Note: This function is essentially a wrapper for the gammatonefir 
%   function in the LTFAT toolbox. The output, H, will be an N-by-M matrix 
%   of filters, where M is the number of elements in FC, rather than the 
%   length M cell array returned by gammatonefir.
%   
%   Needs: LTFAT toolbox.
%
%   See also GAMMATONEFIR, GETPATTERSONFILTERS.

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

narginchk(3,4);

if nargin < 4
    fType = 'complex';
end

validateattributes(fType,{'char'},{'scalartext','nonempty'},...
    'getGammatoneFilters','type of impulse response',4);

switch lower(fType)
    case 'real'
        hcell = gammatonefir(fc,Fs,IRLen,'real');
    case 'complex'
        hcell = gammatonefir(fc,Fs,IRLen);
    otherwise
        error('The fourth input may be ''real'' or ''complex'' only.')
end

numfc = length(fc);
h = zeros(IRLen,numfc);
for ii = 1:numfc
    tempLen = length(hcell{ii}.h);
    h(1:tempLen,ii) = hcell{ii}.h;
end

end
