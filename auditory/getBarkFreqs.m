function [fC,fL,fU,bw] = getBarkFreqs()
%GETBARKFREQS Bark scale frequencies
%   [FC,FL,FU,BW] = GETBARKFREQS returns the 24 center frequencies, FC, of 
%   the Bark scale, along with the band edges, fL and fU, and bandwidths, 
%   BW, all in Hz. Each of the outputs are returned as column vectors.

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

%   Ref:
%       [1]. Fastl and Zwicker (2007) - Psychoacoustics: Facts and Models

narginchk(0,0);

fC = [50, 150, 250, 350, 450, 570, 700, 840, 1000, 1170, 1370, 1600,...
    1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800, 7000, 8500, 10500,...
    13500].';
fL = [0, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720,...
    2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000].';
fU = [100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, ...
    2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, ...
    15500].';
bw = diff([0, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480,...
    1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500,...
    12000, 15500]).';

end
