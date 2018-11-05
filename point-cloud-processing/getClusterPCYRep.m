function [cPCYData,YCoeffs,YMat,degVec] = getClusterPCYRep(cPCData,maxN,...
    varargin)
%GETCLUSTERPCYREP Compute spherical harmonic representation of a clustered
%point cloud.
%   [cPCYData,YCoeffs,YMat,degVec] = GETCLUSTERPCYREP(cPCData,maxN) 
%   computes a real-valued spherical harmonic representation of clustered 
%   point clouds, cPCData. To generate cPCData, see CLUSTERPC. maxN must be 
%   a non-negative scalar corresponding to the maximum spherical harmonic 
%   degree to use. YMat, YCoeffs, and cPCYData are returned as cell arrays
%   with the same number of elements as cPCData. The elements of YMat and 
%   YCoeffs correspond to the spherical harmonics and their corresponding 
%   coefficients, respectively. cPCYData are the clustered point clouds 
%   that are reconstructed from spherical harmonics and has the same format 
%   as cPCData. A vector of spherical harmonic degrees (subject to a
%   maximum of maxN) used to represent each clustered point cloud is
%   returned in degVec.
%
%   ___ = GETCLUSTERPCYREP(...,'Name','Value') allows specification of
%   Name-Value pairs. For options, see GETPCYREP.
%
%   See also CLUSTERPC, GETPCYREP.

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

narginchk(2,6);

% Parse and verify inputs
inputs = parseGETCLUSTERPCYREPInputs(cPCData,maxN,varargin);

% Extract parsed inputs
cPCData = inputs.cPCData;
maxN = inputs.maxN;
typeVal = inputs.TYPE;
csPhaseFlag = inputs.CSPHASE;

% Main computation begins

numCs = length(cPCData);
cPCYData = cell(numCs,1);
YCoeffs = cell(numCs,1);
YMat = cell(numCs,1);
degVec = zeros(numCs,1);
for ii = 1:numCs
    [cPCYData{ii,1},YCoeffs{ii,1},YMat{ii,1},degVec(ii)] = ...
        getPCYRep(cPCData{ii},maxN,'TYPE',typeVal,'CSPHASE',csPhaseFlag,...
        'checkMaxD',1);
end

% Main computation ends

end

function inputs = parseGETCLUSTERPCYREPInputs(cPCData,maxN,opts)
%PARSEGETCLUSTERPCYREPINPUTS Parse and verify inputs for the 
%getClusterPCYRep function.

p = inputParser;

% Required inputs
addRequired(p,'cPCData',@(x)validateattributes(x,{'cell'},{'nonempty'},...
    'getClusterPCYRep','cPCData',1));
addRequired(p,'maxN',@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','integer','nonnegative'},'getClusterPCYRep','maxN',2));

% Optional inputs
addParameter(p,'TYPE','real',@(x)validateattributes(x,{'char'},...
    {'nonempty'},'getClusterPCYRep','TYPE'));
addParameter(p,'CSPHASE',0,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','integer','nonnegative','<=',1},...
    'getClusterPCYRep','CSPHASE'));

p.CaseSensitive = false;
p.FunctionName = 'getClusterPCYRep';

parse(p,cPCData,maxN,opts{:});

inputs = p.Results;

end
