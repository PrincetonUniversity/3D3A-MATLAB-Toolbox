function y = clip(x,xlim,yval)
%CLIP Limit the range of values of a signal.
%   y = CLIP(x) gives x clipped to be between -1 and +1.
%   
%   y = CLIP(x,[xmin,xmax]) gives x for xmin <= x <= xmax, xmin for
%   x < xmin, and xmax for x > xmax.
%   
%   y = CLIP(x,[xmin,xmax],[ymin,ymax]) gives vmin for x < xmin and ymax
%   for x > xmax.
%   
%   y = CLIP(x,xmax) and CLIP(x,xmax,ymax) apply an upper bound only.
%   
%   Note: for xmax > 0, sign(x).*CLIP(abs(x),xmax,ymax) is equivalent to
%   CLIP(x,[-xmax,xmax],[-ymax,ymax]).

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2017 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   ==============================================================================

% Needs at least 1 input argument
narginchk(1,3);

% Gives x clipped to be between -1 and +1 by default
if nargin < 2 || isempty(xlim)
    xlim = [-1 1];
else
    if numel(xlim) > 2
        error('Too many elements of limit vector.');
    end
end

% Clips to limits by default
if nargin < 3 || isempty(yval)
    yval = xlim;
else
    if numel(yval) ~= numel(xlim)
        error('Limit and value vectors must have same number of elements.');
    end
end

y = x;
switch numel(xlim)
    case 1 % Upper bound only
        xmax = xlim;
        ymax = yval;
    case 2 % Upper and lower bounds
        xmin = min(xlim);
        ymin = min(yval);
        xmax = max(xlim);
        ymax = max(yval);
        y(y<xmin) = ymin;
end
y(y>xmax) = ymax;

end