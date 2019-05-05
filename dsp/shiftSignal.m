function out = shiftSignal(x,s)
%SHIFTSIGNAL Shift a signal in time.
%   Y = SHIFTSIGNAL(X,S) shifts a signal, X, by S samples. S must be
%   real-valued (fractional samples are allowed). Negative values of S 
%   correspond to advancement in time (i.e. 'left' shift). X may be a 
%   vector or 2D matrix. If X is a matrix, individual signals must be 
%   specified as the columns of the matrix. 
%       If X is a vector, S must be a scalar.
%       If X is a matrix, S may be either a scalar or vector. The length of
%       the vector must equal the number of columns in X.

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

narginchk(2,2);

% Validate inputs
validateattributes(x,{'double'},{'2d','nonempty','nonnan','finite'},...
    'shiftSignal','X',1)
validateattributes(s,{'double'},{'vector','nonempty','nonnan','finite',...
    'real'},'shiftSignal','S',2)

x = shiftdim(x);
[xLen,numCols] = size(x);

if isscalar(s)
    s = s*ones(numCols,1);
else
    s = shiftdim(s);
    if length(s) ~= numCols
        error('Length of S must equal number of columns in X.')
    end
end

nyqIndx = ceil((xLen+1)/2);
out = zeros(xLen,numCols);
for ii = 1:numCols
    X = fft(x(:,ii));
    
    % Apply time shift property of Fourier transform to shift signal.
    if isreal(x(:,ii))
        Xs = X(1:nyqIndx).*exp(-1i*(2.0*pi*(0:(nyqIndx-1)).'/xLen)*s(ii));
        out(:,ii) = ifft(Xs,xLen,1,'symmetric');
    else
        Xs = X.*exp(-1i*(2.0*pi*(0:(xLen-1)).'/xLen)*s(ii));
        out(:,ii) = ifft(Xs,xLen,1);
    end
end

end
