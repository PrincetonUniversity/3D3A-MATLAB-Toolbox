function out = getERBFreqVec(fL,fU,varargin)
%GETERBFREQVEC Vector of ERB-spaced frequencies.
%   C = GETERBFREQVEC(A,B) returns a frequency vector, C, containing
%   frequencies between A and B (both in Hz) with a uniform spacing of 1
%   ERB. The frequencies in C are also in Hz. A must be <= B. 
%
%   C = GETERBFREQVEC(A,B,'bw',N) returns a frequency vector, C, containing
%   frequencies between A and B (both in Hz) with a uniform spacing of N
%   ERBs. The frequencies in C are also in Hz. A must be <= B. 
%
%   C = GETERBFREQVEC(A,B,'numPts',N) returns a frequency vector, C, of 
%   length N, specifying frequencies between A and B (both in Hz) that are 
%   uniformly-spaced on the ERB scale. The frequencies in C are also in Hz.
%   A must be <= B.

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

narginchk(2,4);

% Check inputs
validateattributes(fL,{'double'},{'scalar','nonempty','nonnan','finite',...
    'real'},'getERBFreqVec','A',1);
validateattributes(fU,{'double'},{'scalar','nonempty','nonnan','finite',...
    'real','>=',fL},'getERBFreqVec','B',2);

if ~isempty(varargin) && length(varargin) ~= 2
    error('Unrecognized input.')
end

bwFlag = true;
indx = find(strcmpi(varargin,'numPts'),1);
if ~isempty(indx)
    n = varargin{indx+1};
    validateattributes(n,{'double'},{'scalar','nonempty','nonnan',...
        'finite','integer','positive'},'getERBFreqVec',...
        'N for input ''numPts''',3);
    bwFlag = false;
else
    indx = find(strcmpi(varargin,'bw'),1);
    if ~isempty(indx)
        bw = varargin{indx+1};
        validateattributes(bw,{'double'},{'scalar','nonempty','nonnan',...
            'finite','real','positive'},'getERBFreqVec',...
            'N for input ''bw''',3);
    else
        bw = 1;
    end
end

fL_erb = f2erb(fL); % Frequency, in ERB units, corresponding to fL.
fU_erb = f2erb(fU);
if bwFlag
    fRange = fU_erb-fL_erb;
    n = floor(fRange/bw);
    r = fRange-(n*bw);
    erbVec = linspace(fL_erb+(r/2),fU_erb-(r/2),n+1);
else
    erbVec = linspace(fL_erb,fU_erb,n);
end

out = erb2f(erbVec);

end
