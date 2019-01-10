function h = fftDeconv(y,x,varargin)
%FFTDECONV Deconvolve two causal signals in the frequency domain.
%   h = FFTDECONV(y,x) deconvolves x out of y in the frequency domain to
%   produce h. The length of x must be less than or equal to that of y. If 
%   y and x are matrices, they must have the same number of columns and the 
%   deconvolution is performed between corresponding pairs of column 
%   vectors. If only y is a matrix, then each column of y is deconvolved by 
%   x (which must be either a vector or a scalar). If x is shorter than y, 
%   then x is zero padded on the right to equal the length of y prior to 
%   deconvolution. Consequently, y and x must both be causal. The length of 
%   h equals the length of y.
%
%   h = FFTDECONV(...,'reg',TYPE) allows specification of the type of
%   regularization to perform when inverting x prior to deconvolving out of
%   y. TYPE must be a cell array and can take any of the optional inputs to
%   the COMPUTEINVERSEFILTER function as elements of the cell array. For
%   example, to use 'gardner1994' with a 'dynRange' of 30, specify TYPE as
%   follows: {'gardner1994',{'dynRange',30}}. By default, no regularization 
%   is performed (i.e. TYPE = {'direct'}). 
%
%   See also FFTCONV, COMPUTEINVERSEFILTER.

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

narginchk(2,3);

% Parse and verify inputs
inputs = parseFFTDECONVInputs(y,x,varargin);

% Extract parsed inputs
y = inputs.y;
x = inputs.x;
TYPE = inputs.TYPE;

% If y and/or x are vectors, force them to column vectors.
y = shiftdim(y);
x = shiftdim(x);

% Check compatibility of input dimensions
[yLen,numColsy] = size(y);
[xLen,numColsx] = size(x);

if xLen > yLen
    error('The length of x cannot exceed that of y.')
else
    x = [x;zeros(yLen-xLen,numColsx)];
end

if numColsx == 1 && numColsy > 1
    x = repmat(x,1,numColsy);
elseif numColsx > 1 && numColsy > 1 && numColsy ~= numColsx
    error(['If y and x are both matrices, they must have the same',...
        ' number of columns'])
elseif numColsy == 1 && numColsx > 1
    error('If y is a vector, x cannot be a matrix.')
end

validateattributes(TYPE{1},{'char'},{'scalartext'},'fftDeconv','TYPE{1}')

if length(TYPE) < 2
    xInv = computeInverseFilter(x,TYPE{1});
else
    xInv = computeInverseFilter(x,TYPE{1},TYPE{2});
end

h = fftConv(y,xInv);

end

function inputs = parseFFTDECONVInputs(y,x,opts)
%PARSEFFTDECONVINPUTS Parse and verify inputs for the fftDeconv function.

p = inputParser;

% Required inputs
addRequired(p,'y',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite','real'},'fftDeconv','y',1));
addRequired(p,'x',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite','real'},'fftDeconv','x',2));

% Optional inputs
addParameter(p,'reg',{'direct'},@(x)validateattributes(x,{'cell'},...
    {'nonempty'},'fftDeconv','TYPE'));

p.CaseSensitive = false;
p.FunctionName = 'fftDeconv';

parse(p,y,x,opts{:});

inputs = p.Results;

end
