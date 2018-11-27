function w = raisedcosinewin(L,r)
%RAISEDCOSINEWIN Raised-cosine window.
%   w = RAISEDCOSINEWIN(L,r) returns an L-point raised-cosine window in the 
%   column vector, w. r must be a vector with two elements, r1 and r2, both
%   of which can take values in the range 0 to 0.5. The resulting window is
%   a rectangular window where the first 100*r1 and last 100*r2 percent of 
%   the samples of the window equal parts of a cosine. If r1 = r2, a Tukey
%   window is returned. If r is not specified, r = [0.5,0.5] is assumed.
%
%   See also TUKEYWIN.

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
    r = [0.5,0.5];
end

validateattributes(L,{'double'},{'scalar','nonempty','nonnan','finite',...
    'integer','positive'},'raisedCosineWin','L',1)
validateattributes(r,{'double'},{'vector','nonempty','nonnan','finite',...
    'real','nonnegative','<=',0.5,'numel',2},'raisedCosineWin','r',2)

wStart = tukeywin(L,2*r(1));
wEnd = tukeywin(L,2*r(2));
halfLen = floor((L+1)/2);
w = [wStart(1:halfLen); wEnd(halfLen+1:end)];

end
