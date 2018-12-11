function minPhaseIR = makeMinPhaseIR(ir,METHOD)
%MAKEMINPHASEIR Compute the minimum-phase version of an impulse response. 
%   minPhaseIR = MAKEMINPHASEIR(ir) computes minPhaseIR, the minimum-phase
%   version of the input, ir, using MATLAB's built-in 'rceps' function. The
%   input, ir, may be a vector or matrix. If ir is a matrix, it is assumed
%   the individual impulse responses are stored as column vectors.
% 
%   ___ = MAKEMINPHASEIR(...,METHOD) optionally specifies the method to
%   use to compute minPhaseIR. The options are:
%       (i) 'rceps' - use MATLAB's built-in 'rceps' function (default).
%       (ii) 'hilb' - use the hilbert transform.
%
%   See also RCEPS.

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

if nargin < 2
    METHOD = 'rceps';
end

ir = shiftdim(ir); % Forces ir to a column if it is a row vector.

switch lower(METHOD)
    case 'rceps'
        minPhaseIR = zeros(size(ir));
        numIRs = size(ir,2);
        for ii = 1:numIRs
            [~,minPhaseIR(:,ii)] = rceps(ir(:,ii));
        end
    case 'hilb'
        absTF = abs(fft(ir));
        argTF = imag(hilbert(log(absTF))); % hilbert operates along columns
        minPhaseTF = absTF.*exp(-1i*argTF);
        minPhaseIR = ifft(minPhaseTF,'symmetric');
end

end
