function y = fftConv(h,x,varargin)
%FFTCONV Convolve two signals in the frequency domain.
%   Y = FFTCONV(H,X) circularly convolves H with X in the frequency domain 
%   to produce Y. 
%       If H and X are vectors of length N, Y will be a length-N column
%       vector.
%       If H and X are matrices of size N-by-M, Y will be an N-by-M matrix
%       with convolution being performed between corresponding pairs of
%       column vectors.
%       If only one of H or X is a matrix of size N-by-M, each column of 
%       the matrix is convolved with the other input (which must be a
%       vector of length N) to produce Y with size N-by-M.
%   This command may also be specified as Y = FFTCONV(H,X,'circ');
%
%   Y = FFTCONV(...,'lin') linearly convolves H with X in the frequency 
%   domain to produce Y. If the signals in H and X have lengths P and Q,
%   respectively, then signals in Y will have length P + Q - 1.
%   Linear convolution is implemented as circular convolution after signals 
%   in H and X have been sufficiently zero-padded on the right. 
%   Consequently, it is assumed that both H and X are causal.
%
%   Y = FFTCONV(...,'lin',TRUNC) optionally specifies whether signals in Y 
%   should be truncated to have the same length as signals in H or X. The 
%   options for TRUNC are provided below:
%   ---------------------------------------------
%       TRUNC       Length of signals in Y
%   =============================================
%         0               P + Q - 1
%         1                   P
%         2                   Q
%   ---------------------------------------------
%   where P and Q denote the length of signals in H and X, respectively.
%
%   See also FFTDECONV.

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

% Parse and verify inputs
inputs = parseFFTConvInputs(h,x,varargin);

% Extract parsed inputs
h = inputs.h;
x = inputs.x;
TYPE = inputs.TYPE;
TRUNC = inputs.TRUNC;

h = shiftdim(h); % If vector, force to column.
x = shiftdim(x);
[hLen,numColsh] = size(h);
[xLen,numColsx] = size(x);
if numColsh == 1 && numColsx > 1
    h = repmat(h,1,numColsx);
elseif numColsx == 1 && numColsh > 1
    x = repmat(x,1,numColsh);
elseif numColsx > 1 && numColsh > 1 && numColsh ~= numColsx
    error(['If H and X are both matrices, they must have the same ',...
        'number of columns.'])
end

if strcmpi(TYPE,'lin')
    fftLen = 2^nextpow2(hLen+xLen-1);
    yFull = ifft(fft(h,fftLen).*fft(x,fftLen),fftLen,1,'symmetric');
    switch TRUNC
        case 0
            y = yFull(1:(hLen+xLen-1),:);
        case 1
            y = yFull(1:hLen,:);
        case 2
            y = yFull(1:xLen,:);
        otherwise
            error('Invalid TRUNC specification for option: ''lin''.')
    end
elseif strcmpi(TYPE,'circ')
    if hLen == xLen
        y = ifft(fft(h).*fft(x),'symmetric');
    else
        error('h and x must have the same length for option: ''circ''.')
    end
else
    error(['Invalid specification for type of convolution. Only ',...
        '''circ'' and ''lin'' are valid.'])
end

end

function inputs = parseFFTConvInputs(h,x,opts)
%PARSEFFTCONVINPUTS Parse and verify inputs for the fftConv function.

p = inputParser;

% Required inputs
addRequired(p,'h',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite'},'fftConv','H',1));
addRequired(p,'x',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite'},'fftConv','X',2));

% Optional inputs
addOptional(p,'TYPE','circ',@(x)validateattributes(x,{'char'},...
    {'scalartext','nonempty'},'fftConv','type of convolution'));
if length(opts) > 1
    if strcmpi(opts{1},'circ')
        warning('Ignoring TRUNC specification for option: ''circ''')
    end
    addOptional(p,'TRUNC',0,@(x)validateattributes(x,{'double'},...
        {'nonempty','nonnan','scalar','finite','integer','nonnegative',...
        '<=',2},'fftConv','TRUNC'));
else
    addOptional(p,'TRUNC',0); % Not used, so no validation needed.
end

p.CaseSensitive = false;
p.FunctionName = 'fftConv';

parse(p,h,x,opts{:});

inputs = p.Results;

end
