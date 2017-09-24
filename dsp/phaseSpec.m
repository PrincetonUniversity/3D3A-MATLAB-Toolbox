function outputSpec = phaseSpec(inputIR,type,tol,DIM)
%PHASESPEC Compute phase spectrum in radians.
%   outputSpec = PHASESPEC(inputIR) returns the principal value phase
%   spectrum given an input impulse response (IR). inputIR may be a matrix 
%   of IRs, with the IRs stored as columns. outputSpec has the same 
%   dimensions as inputIR.
%   
%   outputSpec = PHASESPEC(inputIR,type) additionally allows specification
%   of the type of phase spectrum desired. The options are 'pv' from the
%   principal value phase (default) and 'unwrap' for unwrapped phase. If
%   'unwrap' is specified, a tolerance of pi is used for unwrapping the
%   principal value phase spectrum.
%
%   outputSpec = PHASESPEC(inputIR,'unwrap',tol) additionally allows
%   specification of a custom tolerance for unwrapping the principal value
%   phase spectrum. tol must be specified in radians.
%
%   outputSpec = PHASESPEC(...,DIM) specifies the dimension along which to 
%   perform spectrum computation.
%
%   See also UNWRAP, ANGLE, FFT.

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2017 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   ==============================================================================

narginchk(1,4);

inputIRSize = size(inputIR);
if nargin < 4
    DIM = 1;
    if isvector(inputIR)
        inputIR = shiftdim(inputIR);
    end
end
if nargin < 3
    tol = pi;
end
if nargin < 2
    type = 'pv';
end

if strcmpi(type,'unwrap')
    outputSpec = unwrap(angle(fft(inputIR,[],DIM)),tol,DIM);
else
    outputSpec = angle(fft(inputIR,[],DIM));
end

outputSpec = reshape(outputSpec,inputIRSize);

end

