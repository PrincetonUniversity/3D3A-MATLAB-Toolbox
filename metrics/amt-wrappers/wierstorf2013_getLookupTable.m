function lookupTable = wierstorf2013_getLookupTable(hrtfFile, modelToUse)
%WIERSTORF2013_GETLOOKUPTABLE ITD to azimuth lookup table.
%   TABLE = WIERSTORF2013_GETLOOKUPTABLE(PATH,MODEL) returns the lookup
%   TABLE given HRTFs located in the .sofa file at PATH and using the
%   localization MODEL 'dietz2011' (default) or 'lindemann1986'.
%
%   See also WRAP_WIERSTORF2013, ITD2ANGLE_LOOKUPTABLE.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
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

narginchk(1,2);

if nargin < 2 || isempty(modelToUse)
    modelToUse = 'dietz2011';
end

if exist('SOFAload','file') == 2
    % Load HRTF
    HRTFObj = SOFAload(hrtfFile);
else
    error('SOFAload from SOFA API not found.');
end
if strcmpi(HRTFObj.SourcePosition_Type,'cartesian')
    HRTFObj.SourcePosition = sofaC2sofaS(sourcePosition);
    HRTFObj.SourcePosition_Type = 'spherical';
    HRTFObj.SourcePosition_Units = 'degree, degree, metre';
end
HRTFObj.SourcePosition(:,1) = mod(HRTFObj.SourcePosition(:,1)+180,360)-180;

if exist('itd2angle_lookuptable','file') == 2
    % Construct lookup table
    lookupTable = itd2angle_lookuptable(HRTFObj, HRTFObj.Data.SamplingRate, modelToUse);
else
    error('itd2angle_lookuptable from AMTOOLBOX not found.');
end

end