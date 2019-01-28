function [err,err_avg,err_std] = computeITDError(test,ref,TYPE,DIM)
%COMPUTEITDERROR Compute the error between pairs of ITD values.
%   [E,EM,ES] = COMPUTEITDERROR(A,B) computes the signed difference between
%   ITD values in A and B. The error spectrum (i.e. ITD error as a function 
%   of spatial direction) is returned as E along with EM, the average error 
%   (averaging across all directions). The standard deviation spectrum ES 
%   is also returned. Inputs A and B must be of the same size and store ITD 
%   as row vectors.
%
%   [E,EM,ES] = COMPUTEITDERROR(...,TYPE) additionally specifies the type
%   of error to compute. The options are 'rel' (default) which computes the 
%   signed or relative error between A and B, and 'abs' which computes the 
%   absolute error (i.e. l1 norm).
%
%   [E,EM,ES] = ITDERROR(...,DIM) additionally specifies the type of 
%   averaging to perform. The options are:
%       DIM = 1 - perform averaging across entries (typically different 
%       subjects) for a given spatial direction. Averages over rows.
%       DIM = 2 (default) - perform averaging across spatial directions for 
%       each entry (typically a subject). Averages across columns.

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

narginchk(2,4);

if nargin < 4
    DIM = 2;
end

if nargin < 3
    TYPE = 'rel';
end

switch lower(TYPE)
    case 'rel'
        err = ref-test;
    case 'abs'
        err = abs(ref-test);
end

if (DIM > 2) || (DIM < 1)
    error('Invalid specification for DIM.')
else
    err_avg = mean(err,DIM);
    err_std = std(err,0,DIM);
end

end
