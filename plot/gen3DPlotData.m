function [z,p_u,q_u] = gen3DPlotData(in,pos)
%GEN3DPLOTDATA Format data for generating 3D plots.
%   Z = GEN3DPLOTDATA(X,Y) takes an input vector, X, of data corresponding
%   to positions, Y, specified as coordinate pairs, [p,q], and returns a 
%   matrix, Z, of data where the first dimension of the matrix corresponds
%   to unique values of p and the second dimension to unique values of q.
%   If X is a length N vector, Y must be an N-by-2 matrix.
%
%   [Z,P,Q] = GEN3DPLOTDATA(...) additionally returns the unique values of
%   P and Q corresponding to unique values of the first and second
%   coordinates, respectively, specified in Y.
%
%   See also MESHGRID, SCATTEREDINTERPOLANT.

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

narginchk(2,2);

in = shiftdim(in);
p = round(pos(:,1),4);
q = round(pos(:,2),4);
[p_u,~,~] = unique(p);
[q_u,~,~] = unique(q);
[p_g,q_g] = meshgrid(p_u,q_u);
in_interp = scatteredInterpolant(p,q,in,'nearest','nearest');
z = in_interp(p_g(:),q_g(:));
z = reshape(z,size(p_g));

end
