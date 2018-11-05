function [pcYData,YCoeffs,YMat,degVal] = getPCYRep(pcData,maxN,varargin)
%GETPCYREP Compute spherical harmonic representation of a point cloud.
%   [pcYData,YCoeffs,YMat,degVal] = GETPCYREP(pcData,maxN) computes a real-
%   valued spherical harmonic representation of an input point cloud. The 
%   point cloud, pcData, can be specified as an N-by-3 matrix, where N is 
%   the number of points and where each point is specified in SOFA 
%   cartesian coordinates, or as a point cloud object (see Computer Vision 
%   System Toolbox). maxN must be a non-negative scalar corresponding to 
%   the maximum spherical harmonic degree to use. YMat, YCoeffs, and 
%   pcYData are returned. YMat and YCoeffs correspond to the spherical 
%   harmonics and their corresponding coefficients, respectively. pcYData 
%   is the point cloud that is reconstructed from spherical harmonics and 
%   has the same format as pcData. The spherical harmonic degree (subject 
%   to a maximum of maxN) used to represent the point cloud is returned as 
%   degVal (useful when the 'checkMaxD' option is set to 1). By default,
%   degVal is the same as maxN.
%
%   ___ = GETPCYREP(...,'TYPE','real') - use real-valued spherical 
%   harmonics (default).
%
%   ___ = GETPCYREP(...,'TYPE','complex') - use complex-valued spherical
%   harmonics.
%
%   ___ = GETPCYREP(...,'CSPHASE',0) - ignore the Condon-Shortley phase 
%   term (default).
%
%   ___ = GETPCYREP(...,'CSPHASE',1) - include the Condon-Shortley phase 
%   term.
%
%   ___ = GETPCYREP(...,'checkMaxD',1) - limit the max. spherical harmonic
%   degree used to ensure that the calculation of YCoeffs is an 
%   overdetermined problem. This check is not performed by default.
%
%   See also GETYREP, POINTCLOUD.

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

narginchk(2,8);

% Parse and verify inputs
[inputs,pcFlag] = parseGETPCYREPInputs(pcData,maxN,varargin);

% Extract parsed inputs
pcData = inputs.pcData;
maxN = inputs.maxN;
typeVal = inputs.TYPE;
csPhaseFlag = inputs.CSPHASE;
checkMaxDFlag = inputs.checkMaxD;

% Main computation begins

% 1: Convert point cloud from cartesian to spherical coordinates
sphPCData = sofaC2sofaS(pcData);

% 2: Extract radii of points
pcRadii = sphPCData(:,3);

% 3: Represent radii using spherical harmonics
[yPCRadii,YCoeffs,YMat,degVal] = getYRep(pcRadii.',{'YDirs',pcData,'D',...
    maxN},'TYPE',typeVal,'CSPHASE',csPhaseFlag,'checkMaxD',checkMaxDFlag);

% 4: Reconstruct spherical harmonic representation of point cloud using
% spherical harmonic representation of point radii computed in step 3.
sphYPCData = [sphPCData(:,1:2),yPCRadii.'];

% 5. Convert reconstructed point cloud from spherical to cartesian
% coordinates.
pcYData = sofaS2sofaC(sphYPCData);

% Main computation ends

% Make pcYData have the same format as pcData
if pcFlag
    pcYData = pointCloud(pcYData);
    pcYData.Normal = pcnormals(pcYData);
end

end

function [inputs,pcFlag] = parseGETPCYREPInputs(pcData,maxN,opts)
%PARSEGETPCYREPINPUTS Parse and verify inputs for the getPCYRep function.

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
    'getPCYRep','pcData',1));
addRequired(p,'maxN',@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','integer','nonnegative'},'getPCYRep','maxN',2));

% Optional inputs
addParameter(p,'TYPE','real',@(x)validateattributes(x,{'char'},...
    {'nonempty'},'getPCYRep','TYPE'));
addParameter(p,'CSPHASE',0,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','integer','nonnegative','<=',1},'getPCYRep',...
    'CSPHASE'));
addParameter(p,'checkMaxD',0,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','integer','nonnegative','<=',1},'getPCYRep',...
    'checkMaxD'));

p.CaseSensitive = false;
p.FunctionName = 'getPCYRep';

parse(p,pcData,maxN,opts{:});

inputs = p.Results;

end
