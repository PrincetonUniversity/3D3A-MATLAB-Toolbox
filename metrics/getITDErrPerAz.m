function [mE,sdE,azVec,aE,ppts] = getITDErrPerAz(ITDErr,sPos,BINFLAG)
%GETITDERRPERAZ ITD error per interaural azimuth.
%   [ME,SE,AZ,AE,PP] = GETITDERRPERAZ(E,S) takes a length-N vector of ITD 
%   errors in seconds and an N-by-3 matrix of corresponding source 
%   positions in SOFA cartesian coordinates and returns the mean absolute 
%   ITD error, ME, in seconds at each interaural azimuth, AZ. Also returned 
%   is the standard deviation. The vector of interaural azimuths, AZ, is 
%   computed as azimuths on the horizontal plane. Positions in S that do 
%   not correspond to one of these azimuths are binned with the nearest
%   azimuth (computed in an l2 sense). ME, SE, and AZ are all column
%   vectors of the same length, M. Also returned is an M-by-2 matrix of
%   signed azimuth errors due to binning and a length-M vector PP that 
%   contains the percentage of sample points used to compute ME and SE at 
%   each AZ.
%
%   [ME,SE,AZ,AE,PP] = GETITDERRPERAZ(E,S,false) prevents azimuth binning 
%   as described above.

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

narginchk(2,3);

% Check inputs
validateattributes(ITDErr,{'numeric'},{'vector','real','finite'},...
    'getITDErrPerAz','E',1);
ITDErr = shiftdim(ITDErr); % If vector, force to column vector
numPos = length(ITDErr);
validateattributes(sPos,{'numeric'},{'2d','nonnan','real','finite',...
    'nonempty','size',[numPos,3]},'getITDErrPerAz','S',2);

if nargin < 3
    BINFLAG = true;
end

validateattributes(BINFLAG,{'logical'},{'scalar','nonnan','nonempty'},...
    'getITDErrPerAz','flag for azimuth binning',3);

sDirInt = round(sofaC2cipicI(sPos,'flipAz'),1);
if BINFLAG
    % Get vector of interaural azimuths as azimuths on the horizontal plane
    [~,hpIndxs] = getHPData(sPos.',sPos);   
else
    [~,hpIndxs] = unique(sDirInt(:,[1,2]),'stable','rows');
end
azVec = shiftdim(sDirInt(hpIndxs,1));

% Recompute directions
sDir_b = sDirInt;
azErr = zeros(numPos,1);
for ii = 1:numPos
    sDir_b(ii,1) = findNearest(azVec,sDirInt(ii,1),'l2');
    azErr(ii,1) = sDirInt(ii,1)-sDir_b(ii,1);
end

% Compute mean and standard deviation
numAz = length(azVec);
mE = zeros(numAz,1);
sdE = zeros(numAz,1);
aE = zeros(numAz,2);
ppts = zeros(numAz,1);
for ii = 1:numAz
    azIndxs = find(sDir_b(:,1) == azVec(ii));
    mE(ii,1) = mean(ITDErr(azIndxs));
    sdE(ii,1) = std(ITDErr(azIndxs));
    aE(ii,1) = max(azErr(azIndxs));
    aE(ii,2) = min(azErr(azIndxs));
    ppts(ii,1) = length(azIndxs)*100/numPos;
end

end
