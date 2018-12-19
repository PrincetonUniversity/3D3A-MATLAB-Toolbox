function output = shiftSignal(input,shift)
%SHIFTSIGNAL Shifts a signal in time.
%   output = SHIFTSIGNAL(input,shift) shifts a signal specified by the 
%   sample amount specified in shift. Fractional samples may be specified. 
%   Negative values of shift correspond to advancement in time (i.e. 'left'
%   shift). input may be a vector or 2D matrix. If input is a matrix, the 
%   signals must be specified as the columns of the matrix. shift can be a 
%   scalar or vector (if input is a matrix), with the length of shift equal
%   to the number of columns in input.

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

narginchk(2,2);

input = shiftdim(input);
[signalLen,numCols] = size(input);

if isscalar(shift)
    shift = shift*ones(numCols,1);
else
    shift = shiftdim(shift);
    if length(shift) ~= numCols
        error('Length of shift must equal number of columns in input.')
    end
end

nyqIndx = ceil((signalLen+1)/2);
output = zeros(signalLen,numCols);
for ii = 1:numCols
    inputTF = fft(input(:,ii));
    
    % Apply time shift property of Fourier transform.
    shiftedInput = inputTF(1:nyqIndx).*...
        exp(-1i*(2.0*pi*(0:(nyqIndx-1)).'/signalLen)*shift(ii));
    output(:,ii) = ifft(shiftedInput,signalLen,1,'symmetric');
end

end
