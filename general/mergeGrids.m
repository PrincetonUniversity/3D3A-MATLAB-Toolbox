function [r, a, indx] = mergeGrids(r1, r2, thresh)
%MERGEGRIDS Add missing points to a spherical grid.
%   R = MERGEGRIDS(R1,R2) takes in two grids of positions, R1 and R2, the
%   rows of which correspond to position vectors (given in Cartesian
%   coordinates), and creates a combined grid, R, which contains all points
%   from R1 in addition to those points from R2 that are not already
%   contained in R.
%   
%   R = MERGEGRIDS(R1,R2,T) adds those points from R2 that are at least T
%   degrees away from the closest points in R1. The default threshold is 1
%   degree.
%   
%   [R,A,I] = MERGEGRIDS(...) additionally returns the grid of added
%   positions, A, as well as the corresponding list of row indices of R2,
%   such that R = [R1;A] and A = R2(I).
%
%   See also FINDNEAREST, UNIQUE.

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

if nargin < 3
    thresh = 1; % threshold of 1 degree by default
end

r1Mag = sqrt(dot(r1,r1,2));
r2Mag = sqrt(dot(r2,r2,2));
r2 = r2 * mean(r1Mag) / mean(r2Mag);
u1 = r1 * diag(1./r1Mag);
u2 = r2 * diag(1./r2Mag);

numDirs = size(r2,1);
mask = false(numDirs,1);
for ii = 1:numDirs
    v = u2(ii,:);
    [u,~] = findNearest(u1,v,'angular',1);
    dist = (180/pi) * acos(dot(u,v));
    mask(ii) = (dist > thresh);
end
temp = (1:numDirs).';
indx = temp(mask);

a = r2(indx);
r = [r1; a];

end