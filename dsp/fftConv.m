function y = fftConv(h,x,varargin)
%FFTCONV Convolve two CAUSAL signals in the frequency domain.
%   y = FFTCONV(h,x) convolves h with x in the frequency domain to
%   produce y. By default, circular convolution is performed. The length of
%   y is the length of the longer of h and x. If h and x are matrices, they
%   must have the same number of columns and the convolution is performed 
%   between corresponding pairs of column vectors. If only one of h or x is
%   a matrix, then each column of the matrix is convolved with the other 
%   input (which must not be a matrix). If h and x have different lengths, 
%   the shorter amongst h and x is zero padded on the right to equal the 
%   length of the longer prior to convolution. Consequently, h and x must 
%   both be causal.
%
%   y = FFTCONV(...,TYPE) optionally specifies the type of convolution
%   to perform. The two options are:
%       'lin' - linear convolution
%       'circ' - circular convolution (default)
%   If 'circ' is specified, the shorter amongst h and x is zero padded on
%   the right to equal the length of the longer prior to convolution. The
%   length of y is the length of the longer of h and x.
%   If 'lin' is specified, the length of y is the length of h plus the
%   length of x minus 1.
%
%   y = FFTCONV(...,TRUNC) optionally specifies whether the output, y,
%   should be truncated to have the same length as h or x, if TYPE is 
%   specified as 'lin'. This option does not affect TYPE = 'circ'. The 
%   options for TRUNC are:
%       0 - length of h + length of x - 1 (default)
%       1 - length of h
%       2 - length of x
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
inputs = parseFFTCONVInputs(h,x,varargin);

% Extract parsed inputs
h = inputs.h;
x = inputs.x;
TYPE = inputs.TYPE;
TRUNC = inputs.TRUNC;

h = shiftdim(h);
x = shiftdim(x);
[hLen,numColsh] = size(h);
[xLen,numColsx] = size(x);
if numColsh == 1 && numColsx > 1
    h = repmat(h,1,numColsx);
elseif numColsx == 1 && numColsh > 1
    x = repmat(x,1,numColsh);
elseif numColsx > 1 && numColsh > 1 && numColsh ~= numColsx
    error(['If h and x are both matrices, they must have the same ',...
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
            error('Invalid TRUNC specification')
    end
elseif strcmpi(TYPE,'circ')
    maxLen = max([xLen,hLen]);
    fftLen = 2^nextpow2(maxLen);
    y = ifft(fft(h,fftLen).*fft(x,fftLen),fftLen,1,'symmetric');
    y = y(1:maxLen,:);
else
    error('Invalid TYPE specification')
end

end

function inputs = parseFFTCONVInputs(h,x,opts)
%PARSEFFTCONVINPUTS Parse and verify inputs for the fftConv function.

p = inputParser;

% Required inputs
addRequired(p,'h',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite','real'},'fftConv','h',1));
addRequired(p,'x',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite','real'},'fftConv','x',2));

% Optional inputs
addOptional(p,'TYPE','circ',@(x)validateattributes(x,{'char'},...
    {'scalartext','nonempty'},'fftConv','TYPE'));
addOptional(p,'TRUNC',0,@(x)validateattributes(x,{'double'},...
    {'nonempty','nonnan','scalar','finite','integer','nonnegative',...
    '<=',2},'fftConv','TRUNC'));

p.CaseSensitive = false;
p.FunctionName = 'fftConv';

parse(p,h,x,opts{:});

inputs = p.Results;

end
