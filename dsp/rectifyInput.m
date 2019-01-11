function out = rectifyInput(in,TYPE)
%RECTIFYINPUT Half or full wave rectification.
%   B = RECTIFYINPUT(A) performs half-wave rectification of the signal(s) 
%   in A. If A is a matrix, each column is half-wave rectified. B has the 
%   same dimensions as A. If the input signal has a DC offset (i.e. non-
%   zero mean), this offset is first removed prior to rectification. This 
%   command may also be specified as:
%       B = RECTIFYINPUT(A,'half');
%
%   B = RECTIFYINPUT(A,'full') performs full-wave rectification.

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

narginchk(1,2);

if nargin < 2
    TYPE = 'half';
end

% Check inputs
validateattributes(in,{'double'},{'2d','nonempty','nonnan','finite',...
    'real'},'rectifyInput','A',1);
validateattributes(TYPE,{'char'},{'scalartext','nonempty'},...
    'rectifyInput','type of rectification',2);

% Get input characteristics
flag = false;
if isrow(in)
    in = in.';
    flag = true;
end
sigLen = size(in,1);

% Remove DC offset prior to rectification
in_mean = mean(in);
in_meanMat = repmat(in_mean,sigLen,1);
in_zm = in-in_meanMat; % Zero-mean version of input signal(s).

% Perform rectification
switch lower(TYPE)
    case 'half'
        out = max(in_zm,0);
    case 'full'
        out = abs(in_zm);
    otherwise
        error(['Unrecognized second input. Only ''half'' and ''full''',...
            'are allowed.'])
end

if flag
    out = out.';
end

end
