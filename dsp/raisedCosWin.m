function winVec = raisedCosWin(winLen,r)
%RAISEDCOSWIN Generate a raised-cosine window.
%   winVec = raisedCosWin(winLen) generates a raised-cosine window of
%   length winLen such that the window onset and roll-off each occurs over
%   half winLen. The output, winVec, is a column vector.
%
%   winVec = raisedCosWin(winLen,r) additionally specifies r, which is a
%   1x2 vector where the first element is the fraction of winLen for the 
%   window onset and the second element is the fraction of winLen for the 
%   window roll-off. The default value for r is [0.5,0.5].

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
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

narginchk(1,2);

if nargin == 1
    r = [0.5,0.5];
end

x = (0:1/winLen:(1-(1/winLen)))';
winVec = ones(winLen,1);
r = clip(r,[0,1]);
if r(1) ~= 0
    winVec(1:round(winLen*r(1))) = 0.5*(1+cos((pi()/r(1))*...
        (x(1:round(winLen*r(1)))-r(1))));
end
r(2) = clip(r(2),[0,1-r(1)]);
if r(2) ~= 0
    winVec(winLen-round(winLen*r(2))+1:end) = 0.5*(1+cos((pi()/r(2))*...
       (x(winLen-round(winLen*r(2))+1:end)-1+r(2))));
end

end

