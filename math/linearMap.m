function Y = linearMap(X, XRANGE, YRANGE)
%LINEARMAP Linearly map data from one range to another.
%   Y = LINEARMAP(X,[XMIN,XMAX],[YMIN,YMAX]) maps the data point(s) X from
%   the range [XMIN,XMAX] to the range [YMIN,YMAX] and returns the
%   resulting data point(s) Y.
%
%   If empty or omitted, [YMIN,YMAX] are taken to be [0,1]. Similarly, if
%   empty or omitted, [XMIN,XMAX] are taken to be [MAX(X(:)),MIN(X(:))].
%
%   If either XMIN or YMIN is omitted, it is assumed to be zero, such that
%   LINEARMAP(X,XMAX,YMAX) is the same as LINEARMAP(X,[0,XMAX],[0,YMAX]).

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolkit.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2018 Princeton University
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

if nargin < 2 || isempty(XRANGE)
    XRANGE = [min(X(:)) max(X(:))];
end
if nargin < 3 || isempty(YRANGE)
    YRANGE = [0 1];
end

if numel(XRANGE) == 1
    XRANGE = [0 XRANGE];
end
if numel(YRANGE) == 1
    YRANGE = [0 YRANGE];
end

XMIN = min(XRANGE);
XMAX = max(XRANGE);
YMIN = min(YRANGE);
YMAX = max(YRANGE);

COEFF = (YMAX - YMIN) / (XMAX - XMIN);
Y = COEFF * (X - XMIN) + YMIN;

end