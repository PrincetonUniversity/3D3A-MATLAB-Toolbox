function out = shiftSignal(x,s,varargin)
%SHIFTSIGNAL Shift a signal in time.
%   Y = SHIFTSIGNAL(X,S) shifts a signal, X, by S samples. S must be
%   real-valued (fractional samples are allowed). Negative values of S 
%   correspond to advancement in time (i.e. 'left' shift). If S is not an 
%   integer, a truncated sinc function is used to interpolate the signal 
%   when applying the shift.
%       If X is a vector, S must be a scalar.
%       If X is a matrix, S may be either a scalar or vector. The length of
%       the vector must equal the number of columns in X.
%
%   Y = SHIFTSIGNAL(X,S,METHOD) specifies the method to use to apply the
%   shift when the shift amount is not an integer. METHOD must be specified
%   as a cell array and can take two options:
%       1. {'sinc'} (default) - a truncated sinc function is used to 
%       interpolate the signal when applying the shift.
%       2. {'lagrange',M} - lagrange interpolation of order M is used to
%       interpolate the signal when applying the shift. For example, if M 
%       is 1, linear interpolation is used.
%   If the shift amount, S, is an integer, both methods give the same
%   result.
%
%   Needs: DSP System Toolbox

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
%   Copyright (c) 2019 Princeton University
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

% Check input count
narginchk(2,3);

% Validate required inputs
validateattributes(x,{'double'},{'2d','nonempty','nonnan','finite'},...
    'shiftSignal','X',1)
validateattributes(s,{'double'},{'vector','nonempty','nonnan','finite',...
    'real'},'shiftSignal','S',2)

% Validate optional inputs
if nargin < 3
    METHOD = {'sinc'};
else
    METHOD = varargin{1};
    validateattributes(METHOD,{'cell'},{'2d','nonempty'},'shiftSignal',...
        'METHOD',3)
end

% Define some required variables
x = shiftdim(x);
[xLen,numCols] = size(x);

% Prepare shift vector
if isscalar(s)
    s = s*ones(numCols,1);
else
    s = shiftdim(s);
    if length(s) ~= numCols
        error('Length of S must equal number of columns in X.')
    end
end

% Shift signal
out = zeros(xLen,numCols); % Initialize output
switch lower(METHOD{1})
    case 'sinc'
        nyqIndx = ceil((xLen+1)/2);
        for ii = 1:numCols
            if floor(s(ii)) == s(ii) % If shift amount is an integer
                out(:,ii) = circshift(x(:,ii),s(ii));
            else % If shift amount is not an integer
                X = fft(x(:,ii));
                % Apply time shift property of Fourier transform.
                if isreal(x(:,ii))
                    Xs = X(1:nyqIndx).*exp(-1i*(2.0*pi*(0:(nyqIndx-1)).'...
                        /xLen)*s(ii));
                    out(:,ii) = ifft(Xs,xLen,1,'symmetric');
                else
                    Xs = X.*exp(-1i*(2.0*pi*(0:(xLen-1)).'/xLen)*s(ii));
                    out(:,ii) = ifft(Xs,xLen,1);
                end
            end
        end
    case 'lagrange'
        if numel(METHOD) ~= 2
            error(['For the lagrange interpolation method, order, M, ',...
                'must also be specified.'])
        end
        M = METHOD{2};
        for ii = 1:numCols
            if floor(s(ii)) == s(ii) % If shift amount is an integer
                out(:,ii) = circshift(x(:,ii),s(ii));
            else % If shift amount is not an integer
                s_p1 = floor(s(ii));
                out(:,ii) = circshift(x(:,ii),s_p1);
                s_p2 = s(ii) - s_p1;
                s_sec = s_p2; % Fractional part of shift amount in secs
                d = fdesign.fracdelay(s_sec,'N',M);
                fdFilt = design(d,'lagrange');
                out(:,ii) = filter(fdFilt,out(:,ii));
            end
        end
    otherwise
        error('Invalid METHOD specification.')
end

end
