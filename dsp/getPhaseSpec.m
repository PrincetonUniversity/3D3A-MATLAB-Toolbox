function outputSpec = getPhaseSpec(inputIR,TYPE,TOL,DIM)
%GETPHASESPEC Compute phase spectrum.
%   outputSpec = GETPHASESPEC(inputIR) returns the principal value phase
%   spectrum given an input impulse response (IR). inputIR may be a matrix 
%   of IRs, with the IRs stored as columns. outputSpec has the same 
%   dimensions as inputIR.
%   
%   ___ = GETPHASESPEC(...,TYPE) additionally allows specification of the 
%   type of phase spectrum desired. The options are 'pv' for the principal 
%   value phase (default) and 'unwrap' for unwrapped phase. If 'unwrap' is 
%   specified, a tolerance of pi is used for unwrapping the principal value 
%   phase spectrum.
%
%   ___ = GETPHASESPEC(...,'unwrap',TOL) additionally allows specification 
%   of a custom tolerance for unwrapping the principal value phase 
%   spectrum. tol must be specified in radians. See UNWRAP for more 
%   information.
%
%   ___ = GETPHASESPEC(...,DIM) specifies the dimension along which to 
%   perform spectrum computation.
%
%   See also UNWRAP, ANGLE, GETMAGSPEC.

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

narginchk(1,4);

validateattributes(inputIR,{'double'},{'2d','nonempty','nonnan',...
    'finite','real'},'getPhaseSpec','inputIR',1)

if nargin < 4
    if isvector(inputIR)
        [~,DIM] = max(size(inputIR));
    else
        DIM = 1;
    end
end
validateattributes(DIM,{'double'},{'scalar','nonempty','nonnan',...
    'finite','integer','positive','<=',ndims(inputIR)},'getPhaseSpec',...
    'DIM')

if nargin < 3 || isempty(TOL)
    TOL = pi;
end

if nargin < 2 || isempty(TYPE)
    TYPE = 'pv';
end
validateattributes(TYPE,{'char'},{'scalartext','nonempty'},...
    'getPhaseSpec','TYPE')

switch lower(TYPE)
    case 'pv'
        outputSpec = angle(fft(inputIR,[],DIM));
    case 'unwrap'
        validateattributes(TOL,{'double'},{'scalar','nonempty','nonnan',...
            'finite','real'},'getPhaseSpec','DIM')
        outputSpec = unwrap(angle(fft(inputIR,[],DIM)),TOL,DIM);
    otherwise
        error('Invalid TYPE specification.');
end

end
