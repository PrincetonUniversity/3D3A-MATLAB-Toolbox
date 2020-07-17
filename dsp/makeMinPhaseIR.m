function minPhaseIR = makeMinPhaseIR(ir,varargin)
%MAKEMINPHASEIR Compute the minimum-phase version of an impulse response. 
%   Y = MAKEMINPHASEIR(X) computes Y, the minimum-phase version of the 
%   input, X, using the 'rceps' function in the Signal Processing Toolbox.
%       If X is a vector, Y will be a column vector with the same length as
%       X.
%       If X is a matrix with dimensions N-by-M, individual impulse 
%       responses must be stored as columns. Y will then be a matrix with
%       dimensions N-by-M.
% 
%   Y = MAKEMINPHASEIR(X,METHOD) optionally specifies the method to use
%   when computing Y. The options are:
%       (i) 'rceps' (default) - use the 'rceps' function in the Signal 
%       Processing Toolbox.
%       (ii) 'hilb' - use the hilbert transform (computed using the
%       'hilbert' function in the Signal Processing Toolbox).
%
%   Needs: Signal Processing Toolbox.
%
%   See also RCEPS, HILBERT.

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

narginchk(1,2);

validateattributes(ir,{'numeric'},{'2d','nonempty','nonnan','finite'},...
    'makeMinPhaseIR','X',1)
ir = shiftdim(ir); % Forces ir to a column if it is a row vector.

if nargin < 2
    METHOD = 'rceps';
else
    METHOD = varargin{1};
    validateattributes(METHOD,{'char'},{'scalartext','nonempty'},...
        'makeMinPhaseIR','METHOD',2)
end

minmag = 10*eps; % To deal with ill-conditioned values.
switch lower(METHOD)
    case 'rceps'
        if ~isreal(ir)
            error(['Input X is not real. X must be real when METHOD is',...
                ' ''rceps''.'])
        end
        
        minPhaseIR = zeros(size(ir));
        numIRs = size(ir,2);        
        for ii = 1:numIRs
            [~,minPhaseIR(:,ii)] = rceps(ir(:,ii));
        end
    case 'hilb'
        absTF = getMagSpec(ir);
        absTF(absTF < minmag) = minmag;
        argTF = imag(hilbert(log(absTF))); % hilbert operates along columns
        minPhaseTF = absTF.*exp(-1i*argTF);
        if isreal(ir)
            minPhaseIR = ifft(minPhaseTF,'symmetric');
        else
            minPhaseIR = ifft(minPhaseTF);
        end
    otherwise
        error('Invalid METHOD specification.')
end

end
