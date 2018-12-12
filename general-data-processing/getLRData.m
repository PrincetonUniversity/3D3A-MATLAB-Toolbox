function [lMat,rMat] = getLRData(iMat,cMat,posData,posMat)
%GETLRDATA Convert ipsilateral and contralateral to left and right data. 
%   [lMat,rMat] = GETLRDATA(iMat,cMat,posData,posMat) returns left and  
%   right data matrices (lMat and rMat, respectively), given ipsilateral 
%   and contralateral data matrices (iMat and cMat, respectively) as well 
%   as posData containing position data corresponding to iMat and cMat, and 
%   posMat containing position data corresponding to lMat and rMat. All
%   position data must be in SOFA cartesian coordinates and stored as rows. 
%   Data in iMat and cMat should be stored as column vectors such that the 
%   number of columns in iMat and cMat must each equal the number of rows 
%   in posData and posMat. 
%
%   This function undoes what GETICDATA does.
%
%   See also GETICDATA.

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

narginchk(4,4);

validateattributes(posData,{'double'},{'2d','nonempty','nonnan',...
    'finite','real','size',[NaN,3]},'getLRData','posData',3);
numPos = size(posData,1);
validateattributes(posMat,{'double'},{'2d','nonempty','nonnan',...
    'finite','real','size',[numPos,3]},'getLRData','posMat',4);
validateattributes(iMat,{'double'},{'2d','nonempty','nonnan',...
    'finite','size',[NaN,numPos]},'getLRData','iMat',1);
validateattributes(cMat,{'double'},{'2d','nonempty','nonnan',...
    'finite','size',[NaN,numPos]},'getLRData','cMat',2);

lMat = iMat;
rMat = cMat;
[~,indxs] = ismember(posMat,posData,'rows');
for ii = 1:numPos
    if posMat(ii,2) >= 0
        lMat(:,ii) = iMat(:,indxs(ii));
        rMat(:,ii) = cMat(:,indxs(ii));
    else
        lMat(:,ii) = cMat(:,indxs(ii));
        rMat(:,ii) = iMat(:,indxs(ii));
    end
end

end
