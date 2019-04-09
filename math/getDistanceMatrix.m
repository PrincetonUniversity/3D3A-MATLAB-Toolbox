function D = getDistanceMatrix(X,varargin)
%GETDISTANCEMATRIX Compute matrix of pair-wise distances.
%   D = GETDISTANCEMATRIX(X) computes a matrix of pair-wise distances for
%   the elements in X. If X is an N-by-M matrix, where each row corresponds
%   to the coordinates of a different element in an M-dimensional space,
%   then D will be an N-by-N matrix containing the Euclidean distances 
%   between every permutation of elements in X. Specifically, if x_i, 1 <=
%   i <= N, denote the elements in X, and d_ij, 1 <= i,j <= N denotes the 
%   element in row i and column j of the matrix D, then d_ij is the
%   Euclidean distance between x_i and x_j.
%
%   D = GETDISTANCEMATRIX(X,TYPE) optionally specifies the type of distance
%   to compute. TYPE must be a cell array and can take the following
%   options:
%       1. {'lp',P} - lP norm where P is a positive integer or inf. The 
%       default value of P is 2.
%       2. {'gcd'} - great circle distance.
%
%   See also COMPUTEVECTORNORM, COMPUTEGEODESICMETRIC.

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
inputs = parseGetDistanceMatrixInputs(X,varargin);

% Extract parsed inputs
X = inputs.X;
TYPE = inputs.TYPE;
nR = size(X,1);
D = zeros(nR);

distType = TYPE{1};
validateattributes(distType,{'char'},{'scalartext'},'getDistanceMatrix',...
    'the first entry in TYPE',2);
switch lower(distType)
    case 'lp'
        if length(TYPE) == 1
            P = 2;
        else
            P = TYPE{2};
        end
        validateattributes(P,{'double'},{'scalar','nonempty','nonnan'},...
            'getDistanceMatrix','the value of P for TYPE ''lp''',2);
        
        for ii = 1:nR
            dVec = computeVectorNorm(X-(ones(nR,1)*X(ii,:)),P,2);
            D(ii,:) = dVec.';
        end
    case 'gcd'
        for ii = 1:nR
            gVec = computeGeodesicMetric(X,ones(nR,1)*X(ii,:),2);
            D(ii,:) = gVec.';
        end
    otherwise
        error('Unrecognized first entry in TYPE.')
end

end

function inputs = parseGetDistanceMatrixInputs(X,opts)
%PARSEGETDISTANCEMATRIXINPUTS Parse and verify inputs for the 
%getDistanceMatrix function.

p = inputParser;

% Required inputs
addRequired(p,'X',@(x)validateattributes(X,{'double'},{'2d','nonempty',...
    'nonnan','finite','real'},'getDistanceMatrix','X',1));

% Optional inputs
addOptional(p,'TYPE',{'lp',2},@(x)validateattributes(x,{'cell'},...
    {'vector','nonempty'},'getDistanceMatrix','TYPE',2));

p.CaseSensitive = false;
p.FunctionName = 'getDistanceMatrix';

parse(p,X,opts{:});

inputs = p.Results;

end
