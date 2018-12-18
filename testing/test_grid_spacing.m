% Testing grid spacing

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Joseph G. Tylka <josephgt@princeton.edu>
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


pwOrder = (1:10).';
numPW = (pwOrder + 1).^2;
avgSpacing = zeros(length(pwOrder),1);

for qq = 1:length(pwOrder)
    [r,~] = loadGridFile(['fliege_' int2str(numPW(qq))]);
    
    numDirs = size(r,1);
    distMat = zeros(numDirs);
    for ii = 1:numDirs
        for jj = 1:numDirs
            if ii ~= jj
                distMat(ii,jj) = computeDirectionalError(r(ii,:),r(jj,:));
            end
        end
    end
    
    minAngle = zeros(numDirs,1);
    for ii = 1:numDirs
        distVec = distMat(ii,:);
        minAngle(ii) = min(distVec(distVec~=0));
    end
    
    avgSpacing(qq) = mean(minAngle*180/pi);
end

plot(pwOrder,avgSpacing)
xlabel('plane wave order $\sqrt{Q}-1$','Interpreter','latex')
ylabel('average spacing between adjacent points (\circ)')