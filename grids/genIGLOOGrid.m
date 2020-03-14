function out = genIGLOOGrid(Md)
%GENIGLOOGRID Generate IGLOO spatial sampling grid.
%   B = GENIGLOOGRID() generates an IGLOO spatial sampling grid with
%   resolution 3 following the work of Zhang et al. [1]. The output
%   sampling grid, B, is specified in SOFA cartesian coordinates.
%
%   B = GENIGLOOGRID(A) optionally specifies the resolution, A, to use.

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

%   References:
%       [1] Zhang et al. (2012) - On High-Resolution Head-Related Transfer 
%   Function Measurements: An Efficient Sampling Scheme.

if nargin < 1
    Md = 3;
end

q = (1:(3*2^Md)).';
theta_q = q*(60/(2^Md));
% theta_q = ((2*q)-1)*(30/(2^Md));

Eq = ceil(log2(q))-1;
EMdq = ceil(log2(3*(2^Md)+1-q))-1;

qLen = length(q);
dirCell = cell(qLen,1);
for ii = 1:qLen
    if q(ii) == 1 || q(ii) == (3*(2^Md))
        Vq = 3;
    elseif q(ii) >= 2 && q(ii) <= (2^Md)
        Vq = 9*(2^Eq(ii));
    elseif q(ii) >= (2^Md) + 1 && q(ii) <= 2^(Md + 1)
        Vq = 6*(2^Md);
    elseif q(ii) >= (2^(Md+1)) + 1 && q(ii) <= 3*(2^Md) - 1
        Vq = 9*(2^EMdq(ii));
    end
    
%     if q(ii) == 0 || q(ii) == (3*(2^Md))
%         Vq = 1;
%     elseif q(ii) == 1 || q(ii) == (3*(2^Md) - 1)
%         Vq = 3;
%     elseif q(ii) >= 2 && q(ii) <= (2^Md)
%         Vq = 9*(2^Eq(ii));
%     elseif q(ii) >= ((2^Md) + 1) && q(ii) < 2^(Md + 1)
%         Vq = 6*(2^Md);
%     elseif q(ii) >= (2^(Md+1)) && q(ii) <= (3*(2^Md) - 2)
%         Vq = 9*(2^EMdq(ii));
%     end
    
    v = (0:Vq-1).';
    theta_qVec = repmat(theta_q(ii),Vq,1);
    dirCell{ii,1} = [v*(360/Vq),theta_qVec,ones(Vq,1)];
end

dirMat = cell2mat(dirCell);
dirMat(:,2) = 90-dirMat(:,2);
out = sofaS2sofaC(dirMat);

end
