function Y = logmean(Q,F,FRANGE)
%LOGMEAN Log-scale average of a function.
%   Y = LOGMEAN(Q,F) computes the log-weighted average of Q, whose values
%   are specified as a function of F. F should have uniformly (linearly) 
%   spaced values. Q must have the same number of rows as F.
%
%   Y = LOGMEAN(Q,F,[FMIN,FMAX]) computes the log-weighted average of Q for
%   values of F within FMIN and FMAX.
%
%   See also LOGVAR.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2021 Princeton University
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

% skip DC (F=0) values
Q(F==0,:) = [];
F(F==0) = [];

if nargin==3 && numel(FRANGE)==2
    [~,F1] = findNearest(F,min(FRANGE),'l1');
    [~,F2] = findNearest(F,max(FRANGE),'l1');
    
    if F1==F2
        warning('Only one data point within FRANGE.');
    end
    
    Q = Q(F1:F2,:);
    F = F(F1:F2);
end

% Old weight calculation:
% dF = mean(diff(F));
% W = shiftdim(log((F + dF/2)./(F - dF/2)));

W = shiftdim(1./F);
Y = (W.'*shiftdim(Q))/sum(W);

end
