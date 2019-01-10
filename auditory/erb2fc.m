function fc = erb2fc(erb,n)
%ERB2FC Center frequency corresponding to a specified ERB.
%   B = ERB2FC(A) computes the center frequency, B, in Hz, corresponding to
%   a specified ERB, A, also in Hz. Note that an ERB of 24.7 Hz corresponds
%   to a center frequency of 0.
%
%   B = ERB2FC(A,N) uses the Nth order polynomial approximation given by
%   Moore and Glasberg. Accepts N = 1 (default) or N = 2 only. For N = 2,
%   an ERB of 28.52 Hz corresponds to a center frequency of 0.
%
%   See also FC2ERB, F2ERB, ERB2F.

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
    n = 1;
end

fL = fc2erb(0,n);
switch n
    case 1
        if erb < fL
            warning(['An ERB of %f Hz corresponds to a center',...
                ' frequency of 0 Hz.'],erb)
        end
        
        % The following is based on the formulas on p. 114 in [2].
        % fc = ((erb/24.7)-1)*(1000/4.37);
        
        % The following is based on the Fortran code on p. 135-137 in [2].
        c1 = 24.673;
        c2 = 4.368;
        fc = ((erb/c1)-1)*(1000/c2);
    case 2
        if erb < fL
            warning(['An ERB of %f Hz corresponds to a center',...
                ' frequency of 0 Hz.'],erb)
        end
        fc = 1000*(-93.39 + sqrt(93.39^2 - 4*6.23*(28.52-erb)))/(2*6.23);
    otherwise
        error('No polynomial approximation is known for n = %g',n)
end

end
