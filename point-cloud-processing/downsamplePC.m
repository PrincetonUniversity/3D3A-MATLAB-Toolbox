function outputPC = downsamplePC(inputPC,varargin)
%DOWNSAMPLEPC Downsample a point cloud.
%   outputPC = DOWNSAMPLEPC(inputPC) downsamples an input point cloud so
%   that the resulting point cloud has approx. 5000 points. inputPC can be 
%   an N-by-3 matrix of point cloud data specified in SOFA cartesian 
%   coordinates, or a point cloud object (see Computer Vision System 
%   Toolbox). For more on the algorithm used, see the pcdownsample function
%   in the Computer Vision System Toolbox. outputPC has the same format as
%   inputPC.
%
%   ___ = DOWNSAMPLEPC(...,numPts) optionally specifies the desired 
%   approximate number of points in outputPC.
%
%   See also PCDOWNSAMPLE.

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

narginchk(1,2);

% Parse and verify inputs
[inputs,pcFlag] = parseDOWNSAMPLEPCInputs(inputPC,varargin);

% Extract parsed inputs
inputPC = inputs.pcData;
numPts = inputs.numPts;

% Convert verified inputPC from N-by-3 matrix to point cloud object
inputPC = pointCloud(inputPC);
inputPC.Normal = pcnormals(inputPC);

% Main computation begins

% 1: Approximate initial grid step to use in pcdownsample
YData = getPCYRep(inputPC,0);
approxDia = mean([diff(YData.XLimits),diff(YData.YLimits),...
    diff(YData.ZLimits)]);
approxCurrentGridStep = sqrt(pi*approxDia^2/YData.Count);
initGridStep = approxCurrentGridStep*YData.Count/numPts;
outputPC = pcdownsample(inputPC,'gridAverage',initGridStep);
% 2: Refine grid step approximation
gridStep = initGridStep*sqrt(outputPC.Count/numPts);
outputPC = pcdownsample(inputPC,'gridAverage',gridStep);

% Main computation ends

% Make format of outputPC match that of inputPC
if ~pcFlag
    outputPC = outputPC.Location;
end

end

function [inputs,pcFlag] = parseDOWNSAMPLEPCInputs(pcData,opts)
%PARSEDOWNSAMPLEPCINPUTS Parse and verify inputs for the downsamplePC 
%function.

p = inputParser;

% Required inputs
if strcmpi(class(pcData),'pointCloud')
    pcData = double(pcData.Location);
    pcFlag = 1;
else
    pcFlag = 0;
end
addRequired(p,'pcData',@(x)validateattributes(x,{'pointCloud',...
    'double'},{'2d','nonempty','nonnan','finite','size',[NaN,3]},...
    'downsamplePC','pcData',1));

% Optional inputs
addOptional(p,'numPts',5000,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','positive'},'downsamplePC','numPts'));

p.CaseSensitive = false;
p.FunctionName = 'downsamplePC';

parse(p,pcData,opts{:});

inputs = p.Results;

end
