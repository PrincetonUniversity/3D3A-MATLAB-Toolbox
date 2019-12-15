%% Testing the shiftSignal function

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
%   Copyright (c) 2019 Princeton University
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

%% Specify the signal to be shifted

Fs = 40;
f = 2;
xHalfDur = 1; % Signal duration in seconds
xHalfLen = xHalfDur*Fs;
xLen = round(1.5*xHalfLen);
t = getTimeVec(Fs,xHalfLen);
x = zeros(xLen,1);
x(1:xHalfLen) = sin(2*pi*f*t);

figure();
stem(x);

%% Shift above signal by an integer number of samples

s = 10; % Shift amount in samples
y = shiftSignal(x,s);

figure();
subplot(2,1,1);
stem(x)
subplot(2,1,2);
stem(y)

%% Shift above signal by a fractional number of samples using truncated
% sinc interpolation

s = 10.5; % Shift amount in samples
y = shiftSignal(x,s);

figure();
subplot(2,1,1);
stem(x)
subplot(2,1,2);
stem(y)

%% Shift above signal by a fractional number of samples using cubic 
% lagrange interpolation

s = 10.5; % Shift amount in samples
y = shiftSignal(x,s,{'lagrange',3});

figure();
subplot(2,1,1);
stem(x)
subplot(2,1,2);
stem(y)
