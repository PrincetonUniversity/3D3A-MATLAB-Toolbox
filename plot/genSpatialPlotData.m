function varargout = genSpatialPlotData(dataIn,posIn,varargin)
%GENSPATIALPLOTDATA Format data for generating 2D or 3D spatial plots.
%   Z = GENSPATIALPLOTDATA(D,S) takes a length-N vector, D, of data 
%   corresponding to positions, S, specified as an N-by-2 matrix of 
%   coordinate pairs, and returns a data matrix, Z, with dimensions P-by-Q,
%   where P is the number of unique values in the second column of S and Q
%   is the number of unique values in the first column of S.
%   Nearest-neighbor interpolation and extrapolation of the data in D is
%   performed if necessary. The position values in S are rounded to 5
%   decimal places when extracting unique values to mitigate rounding
%   errors.
%
%   Z = GENSPATIALPLOTDATA(D,S,'posOut',SO) optionally specifies the
%   positions at which the spatial plot data should be returned. A useful
%   function to generate SO might be MAKEPOSMAT.
%
%   Z = GENSPATIALPLOTDATA(D,S,'roundS',false) does not round values in S
%   prior to extracting unique values when 'posOut' and SO are not 
%   specified. If 'posOut' and SO are specified, this command has no effect
%   on the output.
%
%   [Z,X,Y] = GENSPATIALPLOTDATA(__) additionally returns the unique values 
%   in the first column of S as the vector X, and the unique values in the
%   second column of S as the vector Y.
%
%   Needs: MATLAB R2013a or later.
%
%   See also MESHGRID, SCATTEREDINTERPOLANT, UNIQUE.

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
%   Copyright (c) 2019 Princeton University
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

% Check input count
narginchk(2,6);

% Validate required inputs
validateattributes(dataIn,{'double'},{'vector','nonempty','nonnan',...
    'finite','real'},'genSpatialPlotData','D',1);
dataIn = shiftdim(dataIn); % Force dataIn to be a column vector
dataInLen = length(dataIn);
validateattributes(posIn,{'double'},{'2d','nonempty','nonnan','finite',...
    'real','size',[dataInLen,2]},'genSpatialPlotData','S',2);

% Validate optional inputs
indx = find(strcmpi(varargin,'posOut'),1);
if isempty(indx)
    posOutFlag = false;
else
    SO = varargin{indx+1};
    posOutFlag = true;
end

indx = find(strcmpi(varargin,'roundS'),1);
if isempty(indx)
    roundFlag = true;
else
    roundFlag = varargin{indx+1};
end

% Generate interpolation function
if roundFlag
    posInX = round(posIn(:,1),5);
    posInY = round(posIn(:,2),5);
else
    posInX = posIn(:,1);
    posInY = posIn(:,2);
end
interpF = scatteredInterpolant([posInX,posInY],dataIn,'linear','linear');

% Extract unique position values and generate output data matrix
if posOutFlag
    xPos = unique(SO(:,1),'stable');
    yPos = unique(SO(:,2),'stable');
else
    [xPos,~,~] = unique(posInX);
    [yPos,~,~] = unique(posInY);
end
[xPosGrid,yPosGrid] = meshgrid(xPos,yPos);

% Generate output plot data
dataMat = interpF(xPosGrid(:),yPosGrid(:));
dataMat = reshape(dataMat,size(xPosGrid));

% Return outputs
switch nargout
    case 1
        varargout = {dataMat};
    case 3
        varargout = {dataMat,xPos,yPos};
    otherwise
        error('genSpatialPlotData expects request for 1 or 3 outputs.')
end

end
