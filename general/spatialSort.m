function [dataOut,posOut,sI] = spatialSort(dataIn,posIn,varargin)
%SPATIALSORT Sort spatial data based on corresponding position data.
%   Y = SPATIALSORT(X,P) sorts the rows of P such that the elements in its
%   first column are in ascending order, and correspondingly sorts the
%   columns of X. The sorted version of X is returned as the variable Y. X 
%   must have as many columns as rows in P. The sorting of X amounts to a 
%   re-ordering of the columns of X only, since each column in X is assumed
%   to correspond to a row in P. X may be a vector, matrix, or cell array. 
%   P may be a vector or matrix.
%
%   Y = SPATIALSORT(X,P,SO) optionally specifies the sorting order, SO. For 
%   example, SO = [2,1] sorts the rows of P, first in ascending order of
%   the elements in column 2, and then those in column 1. The default value 
%   of SO is 1. For more information, see the COL option in SORTROWS.
%
%   Y = SPATIALSORT(X,P,SO,ROUND) additionally specifies the number of
%   digits to round data in P prior to sorting. The default value is 3. For 
%   more information, see the function ROUND.
%
%   [Y,Q,I] = SPATIALSORT additionally returns, Q, the sorted version of P,
%   and the sort indices, I.
%
%   EXAMPLES OF VALID SYNTAX:
%   
%   Here, A is an N-by-M matrix and B is an M-by-P matrix with P >= 3.
%
%       1. [sA,sB,sI] = spatialSort(A,B); (Assumes DIM = 1)
%       2. [sA,sB,sI] = spatialSort(A,B,[2,1]);
%       3. [sA,sB,sI] = spatialSort(A,B,[2,3,1],2);
%       4. [sA,sB,sI] = spatialSort(A,B,[],1); (Assumes DIM = 1)
%
%   See also SORTROWS.

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

narginchk(2,4);

% Parse and verify inputs
inputs = parseSpatialSortInputs(dataIn,posIn,varargin);

% Extract parsed inputs
dataIn = inputs.dataIn;
posIn = inputs.posIn;
SO = inputs.SO;
ROUND = inputs.ROUND;

if isempty(SO)
    SO = 1;
end

[posOut,sI] = sortrows(round(posIn,ROUND),SO);
dataOut = dataIn(:,sI);

end

function inputs = parseSpatialSortInputs(dataIn,posIn,opts)
%PARSESPATIALSORTINPUTS Parse and verify inputs for the spatialSort
%function.

p = inputParser;

% Required inputs
addRequired(p,'dataIn',@(x)validateattributes(x,{'numeric','cell'},{'2d',...
    'nonempty'},'spatialSort','X',1));
addRequired(p,'posIn',@(x)validateattributes(x,{'numeric'},{'2d',...
    'nonempty','nonnan','finite','size',[size(dataIn,2),NaN]},...
    'spatialSort','P',2));

% Optional inputs
addOptional(p,'SO',1,@(x)validateattributes(x,{'numeric'},{'2d',...
    'integer','positive','<=',size(posIn,2)},'spatialSort','SO',3));
addOptional(p,'ROUND',3,@(x)validateattributes(x,{'numeric'},{'scalar',...
    'nonempty','integer'},'spatialSort','ROUND',4));

p.CaseSensitive = false;
p.FunctionName = 'spatialSort';

parse(p,dataIn,posIn,opts{:});

inputs = p.Results;

end
