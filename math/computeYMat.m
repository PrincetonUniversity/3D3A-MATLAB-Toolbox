function YMat = computeYMat(maxN,dirs,varargin)
%COMPUTEYMAT Compute spherical harmonics.
%   YMat = COMPUTEYMAT(maxN,dirs) returns a matrix of real-valued spherical 
%   harmonics of degrees 0 to maxN for the directions specified in dirs. 
%   maxN must be scalar. dirs must be a p-by-3 matrix of directions 
%   specified in SOFA cartesian coordinates. YMat is then a p-by-(maxN+1)^2 
%   matrix with column 1 corresponding to degree 0, the next 3 columns to 
%   degree 1 and orders -1, 0, and 1, respectively, and so on. The returned 
%   spherical harmonics are orthonormal. The same command may also be
%   specified as:
%       YMat = COMPUTEYMAT(maxN,dirs,'TYPE','real').
%
%   YMat = COMPUTEYMAT(...,'TYPE','complex') returns complex-valued 
%   spherical harmonics instead.
%
%   YMat = COMPUTEYMAT(...,'CSPHASE',0) ignores the Condon-Shortley phase
%   term (default).
%
%   YMat = COMPUTEYMAT(...,'CSPHASE',1) includes the Condon-Shortley phase
%   term.
%
%   See also COMPUTEY.

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

narginchk(2,6);

% Parse and verify inputs
inputs = parseCOMPUTEYMATInputs(maxN,dirs,varargin);

% Extract parsed inputs
maxN = inputs.maxN;
dirs = inputs.dirs;
typeVal = inputs.TYPE;
csPhaseFlag = inputs.CSPHASE;

if strcmpi(typeVal,'real') || strcmpi(typeVal,'complex')
    % Main computation begins
    
    n = (0:maxN)';
    Y = computeY(n,dirs,'TYPE',typeVal,'CSPHASE',csPhaseFlag);
    YMat = cell2mat(Y);
    
    % Main computation ends
else
    error('Unrecognized input. TYPE must be ''%s'' or ''%s.''',...
        'real','complex')
end

end

function inputs = parseCOMPUTEYMATInputs(maxN,dirs,opts)
%PARSECOMPUTEYMATINPUTS Parse and verify inputs for the computeYMat
%function.

p = inputParser;

% Required inputs
addRequired(p,'maxN',@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','nonnan','finite','nonnegative'},'computeYMat','maxN',1));
addRequired(p,'dirs',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite','size',[NaN,3]},'computeYMat','dirs',2));

% Optional inputs

addParameter(p,'TYPE','real',@(x)validateattributes(x,{'char'},...
    {'nonempty'},'computeYMat','TYPE'));
addParameter(p,'CSPHASE',0,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','integer','nonnegative','<=',1},'computeYMat',...
    'CSPHASE'));

p.CaseSensitive = false;
p.FunctionName = 'computeYMat';

parse(p,maxN,dirs,opts{:});

inputs = p.Results;

end
