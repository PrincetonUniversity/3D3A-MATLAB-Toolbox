function outputSpec = getPhaseSpec(inputIR,varargin)
%GETPHASESPEC Compute phase spectrum.
%   S = GETPHASESPEC(X) returns the principal value phase spectrum, S, 
%   given an input impulse response (IR), X. If X is a matrix of IRs, the 
%   IRs must be stored as columns. S will have the same dimensions as X.
%   All phase values are specified in radians.
%   
%   S = GETPHASESPEC(X,TYPE) additionally allows specification of the type
%   of phase spectrum that is computed. The options are 'pv' for the 
%   principal value phase (default) and 'unwrap' for unwrapped phase. If 
%   'unwrap' is specified, a tolerance of pi is used for unwrapping the 
%   principal value phase spectrum.
%
%   S = GETPHASESPEC(X,'unwrap',TOL) additionally allows specification of a
%   custom tolerance for unwrapping the principal value phase spectrum. TOL 
%   must be specified in radians. See the 'unwrap' function for more 
%   information.
%
%   S = GETPHASESPEC(...,DIM) specifies the dimension along which to
%   compute S. TYPE and/or TOL must each be specified as [] to use their
%   default values when specifying a custom value for DIM.
%
%   See also UNWRAP, ANGLE, GETMAGSPEC, GETMAGSPECDB.

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
%   Copyright (c) 2020 Princeton University
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

narginchk(1,4);

validateattributes(inputIR,{'numeric'},{'2d','nonempty','nonnan',...
    'finite'},'getPhaseSpec','X',1)

if nargin < 4
    if isvector(inputIR)
        [~,DIM] = max(size(inputIR));
    else
        DIM = 1;
    end
else
    DIM = varargin{3};
    validateattributes(DIM,{'double'},{'scalar','nonempty','nonnan',...
        'finite','integer','positive','<=',ndims(inputIR)},...
        'getPhaseSpec','DIM',4)
end

if nargin < 3
    TOL = pi;
else
    TOL = varargin{2};
    if isempty(TOL)
        TOL = pi;
    end
end

if nargin < 2
    TYPE = 'pv';
else
    TYPE = varargin{1};
    if isempty(TYPE)
        TYPE = 'pv';
    else
        validateattributes(TYPE,{'char'},{'scalartext','nonempty'},...
            'getPhaseSpec','TYPE',2)
    end
end

switch lower(TYPE)
    case 'pv'
        outputSpec = angle(fft(inputIR,[],DIM));
    case 'unwrap'
        validateattributes(TOL,{'numeric'},{'scalar','nonnan','finite',...
            'real'},'getPhaseSpec','TOL',3)
        outputSpec = unwrap(angle(fft(inputIR,[],DIM)),TOL,DIM);
    otherwise
        error('Invalid TYPE specification.');
end

end
