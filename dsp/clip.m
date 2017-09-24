function y = clip(x,lim,val)
%CLIP Limit the range of values of x.
%   y = CLIP(x) gives x clipped to be between -1 and +1.
%   
%   y = CLIP(x,[min,max]) gives x for min <= x <= max, min for x < min, and
%   max for x > max.
%   
%   y = CLIP(x,[min,max],[vmin,vmax]) gives vmin for x < min and vmax for 
%   x > max.
%   
%   y = CLIP(x,max) and CLIP(x,max,vmax) apply an upper bound only.
%   
%   Note: for max > 0, sign(x).*CLIP(abs(x),max,vmax) is equivalent to
%   CLIP(x,[-max,max],[-vmax,vmax]).

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
if nargin < 1
    error('Not enough input arguments.');
end

% Gives x clipped to be between -1 and +1 by default
if nargin < 2
    lim = [-1 1];
else
    if numel(lim) > 2
        error('Too many elements of limit vector.');
    end
end

% Clips to limits by default
if nargin < 3
    val = lim;
else
    if numel(val) ~= numel(lim)
        error('Limit and value vectors must have same number of elements.');
    end
end

y = x;
switch numel(lim)
    case 1 % Upper bound only
        max = lim;
        vmax = val;
    case 2 % Upper and lower bounds
        min = lim(1);
        vmin = val(1);
        max = lim(2);
        vmax = val(2);
        if min > max
            error('Lower bound must not be greater than upper bound.');
        end
        y(y<min) = vmin;
end
y(y>max) = vmax;

end