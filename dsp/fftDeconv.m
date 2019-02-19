function h = fftDeconv(y,x,varargin)
%FFTDECONV Deconvolve two signals in the frequency domain.
%   H = FFTDECONV(Y,X) circularly deconvolves X out of Y in the frequency 
%   domain to produce H.
%       If X and Y are vectors, they must have lengths M and N,
%       respectively, where M <= N. H will be a column vector of length N.
%       If X and Y are both matrices, they must have dimensions M-by-P and
%       N-by-P, respectively, where M <= N. H will have dimensions N-by-P.
%       If Y is a matrix and X is a vector, each column of Y is deconvolved
%       by X. H will have the same dimensions as Y.
%   If the lengths of the signals in X are shorter than those in Y, then
%   the signals in X are zero padded on the right to equal the lengths of
%   the signals in Y prior to deconvolution. This command may also be
%   specified as: H = FFTDECONV(Y,X,'circ') or H = FFTDECONV(Y,X,[]).
%
%   H = FFTDECONV(Y,X,'padcirc') performs padded circular deconvolution. 
%   This is the same as circular deconvolution as described above, except 
%   signals in both X and Y are zero-padded on the right to have lengths 
%   that are the first power of 2 greater than M+N-1, where M and N are the 
%   lengths of the signals in X and Y, respectively. The signals in H will
%   have the same lengths as those in Y prior to zero padding.
%
%   H = FFTDECONV(Y,X,'lin') performs linear deconvolution such that X,
%   when linearly convolved with H, will produce Y. If the lengths of the 
%   signals in X are M and those in Y are N, where M <= N, then the length
%   of the signals in H will be N-M+1.
%
%   H = FFTDECONV(...,'reg',TYPE) allows specification of the type of
%   regularization to perform when inverting X prior to deconvolving out of
%   Y. TYPE must be a cell array and can take any of the optional inputs to
%   the COMPUTEINVERSEFILTER function as elements of the cell array. For
%   example, to use 'gardner1994' with a 'dynRange' of 30, specify TYPE as
%   follows: {'gardner1994',{'dynRange',30}}. By default, no regularization 
%   is performed (i.e. TYPE = {'direct'}). If 'circ','padcirc' or 'lin' is
%   not specified as a third input, specify the third input as [] before
%   specifying 'reg' and TYPE.
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

narginchk(2,5);

% Parse and verify inputs
inputs = parseFFTDeconvInputs(y,x,varargin);

% Extract parsed inputs
y = inputs.y;
x = inputs.x;
OPT = inputs.OPT;
TYPE = inputs.reg;
if isempty(OPT)
    OPT = 'circ';
end

% If y and/or x are vectors, force them to column vectors.
y = shiftdim(y);
x = shiftdim(x);

% Check compatibility of input dimensions
[yLen,numColsy] = size(y);
[xLen,numColsx] = size(x);

if xLen > yLen
    error('The length of x cannot exceed that of y.')
end

if numColsx == 1 && numColsy > 1
    x = repmat(x,1,numColsy);
elseif numColsx > 1 && numColsy > 1 && numColsy ~= numColsx
    error(['If y and x are both matrices, they must have the same',...
        ' number of columns'])
elseif numColsy == 1 && numColsx > 1
    error('If y is a vector, x cannot be a matrix.')
end

validateattributes(TYPE{1},{'char'},{'scalartext'},'fftDeconv',...
    'first element of ''TYPE'' specification for option ''reg''.')

switch lower(OPT)
    case 'circ'
        x = [x;zeros(yLen-xLen,numColsx)];
        hLen = yLen;
    case 'padcirc'
        padLen = 2^nextpow2(yLen+xLen);
        y = [y;zeros(padLen-yLen,numColsy)];
        x = [x;zeros(padLen-xLen,numColsx)];
        hLen = yLen;
    case 'lin'
        x = [x;zeros(yLen-xLen,numColsx)];
        hLen = yLen-xLen+1;
    otherwise
        error('%s is an unrecognized input.',OPT)
end

if length(TYPE) < 2 % Unregularized inverse
    xInv = computeInverseFilter(x,TYPE{1});
else % Regularized inverse
    xInv = computeInverseFilter(x,TYPE{1},TYPE{2});
end

h = fftConv(y,xInv);
h = h(1:hLen,:);

end

function inputs = parseFFTDeconvInputs(y,x,opts)
%PARSEFFTDECONVINPUTS Parse and verify inputs for the fftDeconv function.

p = inputParser;

% Required inputs
addRequired(p,'y',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite'},'fftDeconv','y',1));
addRequired(p,'x',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite'},'fftDeconv','x',2));

% Optional inputs
addOptional(p,'OPT','circ',@(x)validateattributes(x,{'char','double'},...
    {},'fftDeconv','type of deconvolution',3));
addParameter(p,'reg',{'direct'},@(x)validateattributes(x,{'cell'},...
    {'nonempty'},'fftDeconv','TYPE'));

p.CaseSensitive = false;
p.FunctionName = 'fftDeconv';

parse(p,y,x,opts{:});

inputs = p.Results;

end
