function [T,R,L] = computePCTMatrix(P,Q)
%COMPUTEPCTMATRIX Compute point cloud transformation matrix.
%   T = COMPUTEPCTMATRIX(P,Q) computes an approximate transformation
%   matrix, T, that describes the translation and rotation applied to point
%   cloud P to get to point cloud Q.
%
%   [T,R] = COMPUTEPCTMATRIX(P,Q) additionally returns rotation values as
%   a 3-element row vector [alpha,beta,gamma] containing yaw (alpha),
%   pitch (beta), and roll (gamma) values extracted from the T matrix.
%
%   [T,R,L] = COMPUTEPCTMATRIX(P,Q) additionally returns translation values 
%   as a 3-element row vector [x,y,z] extracted from the T matrix.

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

validateattributes(P,{'double'},{'2d','nonempty','nonnan','finite',...
    'real','size',[NaN,3]},'computePCTMatrix','P',1);
numPts = size(P,1);
validateattributes(Q,{'double'},{'2d','nonempty','nonnan','finite',...
    'real','size',[numPts,3]},'computePCTMatrix','Q',2);

% Express point clouds in homogeneous coordinates

P_h = [P,ones(numPts,1)];
Q_h = [Q,ones(numPts,1)];

% Compute optimal (in least-squares sense) T

T = pinv(P_h)*Q_h;

% Estimate translation and rotation values

if nargout > 1
    alpha = atan2d(T(1,2),T(1,1)); 
    beta = atan2d(-T(1,3),sqrt((T(1,2))^2+(T(1,1))^2)); 
    gamma = atan2d(T(2,3),T(3,3)); 
    R = [gamma,beta,alpha];
end

if nargout > 2
    L = T(4,1:3);
end

end
