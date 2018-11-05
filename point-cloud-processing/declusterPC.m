function pcData = declusterPC(cPCData)
%DECLUSTERPC Generate a single point cloud from a cluster of point clouds.
%   pcData = DECLUSTERPC(cPCData) takes a cell array, cPCData, of clustered
%   point clouds and returns a single point cloud, pcData. pcData is
%   returned as an N-by-3 matrix, where N is the number of points and where 
%   each point is specified in SOFA cartesian coordinates. Each element of
%   cPCData may be specified as an N-by-3 matrix (analogous to pcData), or
%   as a point cloud object (see Computer Vision System Toolbox).
%
%   See also CLUSTERPC.

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

narginchk(1,1);

% Verify input
validateattributes(cPCData,{'cell'},{'nonempty'},'declusterPC','cPCData',1)

% Main computation begins

numCs = length(cPCData); % Compute number of clustered point clouds
for ii = 1:numCs
    if strcmpi(class(cPCData{ii}),'pointCloud')
        cPCData{ii} = double(cPCData{ii}.Location);
    else
        cPCData{ii} = double(cPCData{ii});
    end
end
pcData = cell2mat(shiftdim(cPCData));

% Main computation ends

end
