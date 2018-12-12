function [hPData,hPIndxs] = getHPData(inputData,posMat)
%GETHPDATA Get horizontal plane data.
%   [hPData,hPIndxs] = GETHPDATA(inputData,posMat) returns horizontal
%   plane data given inputData with corresponding positions in posMat.
%   inputData must be a matrix with the data stored as column vectors.
%   posMat must be in cartesian coordinates with the coordinates stored as
%   row vectors. The number of rows in posMat must equal the number of
%   columns in inputData. The indices of posMat corresponding to the
%   horizontal plane are also returned as hPIndxs.
%
%   See also GETLPDATA, GETMPDATA, GETUHDATA.

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

validateattributes(posMat,{'double'},{'2d','nonempty','nonnan',...
    'finite','real','size',[NaN,3]},'getHPData','posMat',2);
numPos = size(posMat,1);
validateattributes(inputData,{'double'},{'2d','nonempty','nonnan',...
    'finite','size',[NaN,numPos]},'getHPData','inputData',1);

% Check to see where the z-coordinate is zero. round is used to prevent
% numerical precision errors from causing z = 0 positions to be missed.  
hPIndxs = round(posMat(:,3),3) == 0; 
hPData = inputData(:,hPIndxs);

end
