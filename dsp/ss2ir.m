function [h,t,x] = ss2ir(A,B,C,D,L,FS)
%SS2IR Finite impulse response from state-space model
%   H = SS2IR(A,B,C,D,L,FS) returns a finite impulse response (FIR) of 
%   length L at sampling rate FS given the state-space matrices A,B,C, and 
%   D. L must be specified in samples and FS in Hz.
%
%   [H,T,X] = SS2IR(A,B,C,D,L,FS) optionally returns the time vector, T,
%   and state-trajectories, X. For more information, see LSIM.
%
%   Needs: Control System Toolbox.
%
%   See also IR2SS, LSIM.

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

narginchk(6,6);

m = size(B,2); % Extract number of inputs
x0 = zeros(size(A,1),1); % Define vector of initial state values
imp = repmat([1;zeros(L-1,1)],1,m); % Generate input signals
t = getTimeVec(FS,L); % Generate time vector (in seconds)

% Simulate impulse response
sys = ss(A,B,C,D,1/FS);
[h_tmp,t,x] = lsim(sys,imp,t,x0);
if D == 0
    h = circshift(h_tmp,-1,1);
else
    h = h_tmp;
end

end
