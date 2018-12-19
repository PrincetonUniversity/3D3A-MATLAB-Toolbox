function [dataOut,posOut,sI] = spatialSort(dataIn,posIn,varargin)
%SPATIALSORT Sort spatial data based on corresponding position data.
%   [dataOut,posOut,sI] = SPATIALSORT(dataIn,posIn) sorts the rows of posIn 
%   such that the elements in its first column are in ascending order, and 
%   correspondingly sorts the columns of dataIn. Sorted versions of dataIn
%   and posIn are returned in dataOut and posOut, respectively. dataIn must 
%   have as many columns as rows in posIn. Sort indices are returned in sI. 
%   The sorting of dataIn amounts to a re-ordering of the columns of dataIn 
%   only, since each column in dataIn corresponds to a row in posIn.
%
%   ___ = SPATIALSORT(...,DIM) optionally specifies the sorting order. For 
%   example, DIM = [2,1] sorts the rows of posIn, first in ascending order
%   of the elements in column 2, and then those in column 1. The default 
%   value of DIM is 1. For more information, see the COL option in the 
%   function SORTROWS.
%
%   ___ = SPATIALSORT(...,ROUND) additionally specifies the number of 
%   digits to round data in posIn prior to sorting. The default value is 
%   3. For more information, see the function ROUND.
%
%   EXAMPLES OF VALID SYNTAX:
%   
%   Here, A is an N-by-M matrix and B is a M-by-P matrix with P >= 3.
%
%       1. [sA,sB,sI] = spatialSort(A,B); (Assumes DIM = 1)
%       2. [sA,sB,sI] = spatialSort(A,B,[2,1]);
%       3. [sA,sB,sI] = spatialSort(A,B,[2,3,1],2);
%       4. [sA,sB,sI] = spatialSort(A,B,[],1); (Assumes DIM = 1)
%
%   See also SORTROWS, ROUND.

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
inputs = parseSPATIALSORTInputs(dataIn,posIn,varargin);

% Extract parsed inputs
dataIn = inputs.dataIn;
posIn = inputs.posIn;
DIM = inputs.DIM;
ROUND = inputs.ROUND;

if isempty(DIM)
    DIM = 1;
end

% Main computation begins
[posOut,sI] = sortrows(round(posIn,ROUND),DIM);
dataOut = dataIn(:,sI);
% Main computation ends

end

function inputs = parseSPATIALSORTInputs(dataIn,posIn,opts)
%PARSESPATIALSORTINPUTS Parse and verify inputs for the spatialSort
%function.

p = inputParser;

% Required inputs
addRequired(p,'dataIn',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite'},'spatialSort','dataIn',1));
addRequired(p,'posIn',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite','size',[size(dataIn,2),NaN]},...
    'spatialSort','posIn',2));

% Optional inputs
addOptional(p,'DIM',1,@(x)validateattributes(x,{'double'},{'2d',...
    'integer','positive','<=',size(posIn,2)},'spatialSort','DIM'));
addOptional(p,'ROUND',3,@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','integer'},'spatialSort','ROUND'));

p.CaseSensitive = false;
p.FunctionName = 'spatialSort';

parse(p,dataIn,posIn,opts{:});

inputs = p.Results;

end
