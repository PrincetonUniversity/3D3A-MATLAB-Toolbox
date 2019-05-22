function posMat = makePosMat(uPosVecs)
%MAKEPOSMAT Construct a matrix of position coordinate groups.
%   GP = MAKEPOSMAT(UP) takes a cell array of vectors of unique coordinates 
%   and exports a matrix of grouped (pairs, triples, etc.) coordinates. UP
%   must be a cell array of length N, with each element being either a
%   scalar or vector. If the elements in UP have lengths M_i, where 
%   i = 1,2,...,N, then GP will be a Q-by-N matrix, where Q = M_1*M_2*...
%   M_(N-1)*M_N.
%
%   EXAMPLE:
%
%       a = 1:3; b = [4;5]; c = 1;
%       d{1} = a; d{2} = b; d{3} = c;
%       posMat = makePosMat(d);

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

narginchk(1,1);

% Validate input
validateattributes(uPosVecs,{'cell'},{'nonempty'},'makePosMat','UP',1)

fprintf('Generating position matrix...');

N = length(uPosVecs); % Count number of unique position vectors
if N == 1
    posMat = shiftdim(uPosVecs{1});
else
    % 1: Compute length of each unique vector in uniquePosVecs
    Q = shiftdim(cellfun(@length,uPosVecs));
    % 2: Compute product of vector lengths to set size of posMat
    P = prod(Q);
    posMat = zeros(P,N);
    % 3: Populate posMat
    cumProdVec = cumprod(flipud(Q)); % Cumulative product of lengths
    repCount = 1; % Required # of repetitions of each element in currentVec
    for ii = 1:N
        partitionVec = zeros(cumProdVec(ii),1);
        currentVec = uPosVecs{N-ii+1};
        
        % Check if each element in uniquePosVecs is a scalar or vector
        validateattributes(currentVec,{'double'},{'vector','nonempty',...
            'finite','nonnan'},'makePosMat',['UP{',num2str(ii),'}'],1)
        
        indx = 1;
        for jj = 1:repCount:cumProdVec(ii)
            hIndx = jj+repCount-1;
            partitionVec(jj:hIndx) = repmat(currentVec(indx),repCount,1);
            indx = indx + 1;
        end
        repCount = cumProdVec(ii); % Update repCount
        numPartitions = P/cumProdVec(ii);
        posMat(:,(N-ii+1)) = repmat(partitionVec,numPartitions,1);
    end
end

fprintf('done!\n');

end
