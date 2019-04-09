function g = computeGeodesicMetric(a,b,DIM)
%COMPUTEGEODESICMETRIC Geodesic metric on S2.
%   G = COMPUTEGEODESICMETRIC(A,B) computes the geodesic metric between two
%   points on the unit sphere, S2, given by the unit vectors, A and B. The
%   unit vectors must be specified in cartesian coordinates. G is specified 
%   in radians in the range [0,pi]. If A and B are matrices specifying
%   multiple points, they must have the same size, and the metric is 
%   computed by assuming that each column corresponds to a different point 
%   (i.e. the components of a position vector are specified across rows).
%
%   G = COMPUTEGEODESICMETRIC(A,B,DIM) specifies the dimension along which
%   the components of the position vector for each point are specified.
%
%   See also DOT.

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

narginchk(2,3);

if nargin < 3
    DIM = find(size(a) > 1,1);
end

g = acos(dot(a,b,DIM));

end
