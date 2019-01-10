function f = erb2f(erb)
%ERB2F Frequency, in Hz, corresponding to input frequency in ERBs.
%   B = ERB2F(A) returns the frequency, B, in Hz, corresponding to the
%   frequency, A, specified in ERB units.
%
%   See also F2ERB, FC2ERB.

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

%   Ref:
%       [1]. Glasberg and Moore (1990) - Derivation of auditory filter 
%       shapes from notched-noise data.

narginchk(1,1);

% The following is from the Fortran code on p. 138 in [1].
c1 = 24.673;
c2 = 4.368;
c3 = 2302.6/(c1*c2);
f = 1000*((10.^(erb/c3))-1)/c2;

% The following is from the formulas on p. 114 in [1].
% f = (exp(erb*24.7*4.37/1000)-1)*1000/4.37;
% f = (10.^(erb/21.4)-1)*(1000/4.37);

end
