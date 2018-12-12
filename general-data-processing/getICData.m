function [iMat,cMat,posData] = getICData(lMat,rMat,posMat)
%GETICDATA Get ipsilateral and contralateral data.
%   [iMat,cMat,posData] = GETICDATA(lMat,rMat,posMat) returns ipsilateral 
%   and contralateral data matrices (iMat and cMat, respectively), given 
%   left-ear and right-ear data matrices (lMat and rMat, respectively) and 
%   position data, posMat. posMat must be in SOFA cartesian coordinates and 
%   stored as rows. Data in lMat and rMat should be stored as column 
%   vectors such that the number of columns in lMat and rMat must each
%   equal the number of rows in posMat. Position data corresponding to the 
%   data in iMat and cMat is also returned as posData.
%
%   See also GETLRDATA.

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

narginchk(3,3);

validateattributes(posMat,{'double'},{'2d','nonempty','nonnan',...
    'finite','real','size',[NaN,3]},'getICData','posMat',3);
numPos = size(posMat,1);
validateattributes(lMat,{'double'},{'2d','nonempty','nonnan',...
    'finite','size',[NaN,numPos]},'getICData','lMat',1);
validateattributes(rMat,{'double'},{'2d','nonempty','nonnan',...
    'finite','size',[NaN,numPos]},'getICData','rMat',2);

lPosIndxs = round(posMat(:,2),3) >= 0;
rPosIndxs = round(posMat(:,2),3) < 0;

lPosMat = posMat(lPosIndxs,:);
rPosMat = posMat(rPosIndxs,:);
posData = [lPosMat;rPosMat];

iMat = [lMat(:,lPosIndxs),rMat(:,rPosIndxs)];
cMat = [rMat(:,lPosIndxs),lMat(:,rPosIndxs)];

end
