function A = computeRBFInterpMatrix(X,varargin)
%COMPUTERBFINTERPMATRIX Compute radial basis function interpolation matrix.
%   A = COMPUTERBFINTERPMATRIX(X) returns an interpolation matrix, A, for
%   distinct input data points, X, using a Gaussian radial basis function
%   (RBF). If X has dimensions N-by-M, then A will have dimensions N-by-N.
%   X must be specified in standard cartesian coordinates, with each row
%   corresponding to a different data point. M therefore corresponds to the
%   dimension of the input space.
%
%   A = COMPUTERBFINTERPMATRIX(X,TYPE) optionally specifies the TYPE of RBF 
%   to use. TYPE must be specified as a cell array and can take the 
%   following options:
%       1. {'Gaussian',C} - Gaussian RBF with C = 20 (default). A different
%       C may be specified, subject to 0 < C < inf. A larger C corresponds 
%       to a narrower 'width' of the Gaussian.
%       2. {'Inverse Multiquadric,C} - Inverse multiquadric RBF with C = 1
%       (default). A different C may be specified, subject to 0 < C < inf.
%       3. {'Multiquadric,C} - Multiquadric RBF with C = 1 (default). A 
%       different C may be specified, subject to 0 < C < inf.
%       4. {'Linear'} - Linear RBF.
%       5. {'TPS'} - Thin plate spline RBF.
%       6. {'Cubic'} - Cubic RBF.

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
inputs = parseComputeRBFInterpMatrixInputs(X,varargin);

% Extract parsed inputs
X = inputs.X;
TYPE = inputs.TYPE;

% Compute distance matrix
dMat = getDistanceMatrix(X);

rbfType = TYPE{1};
validateattributes(rbfType,{'char'},{'scalartext'},...
    'computeRBFInterpMatrix','the first entry in TYPE',2);
switch lower(rbfType)
    case 'gaussian'
        if length(TYPE) == 1
            C = 20;
        else
            C = TYPE{2};
        end
        
        validateattributes(C,{'double'},{'scalar','nonempty','nonnan',...
            'finite','real','positive'},'computeRBFInterpMatrix',...
            'the value of C when TYPE is ''Gaussian''',2);
        A = exp(-C*(dMat.^2));
    case 'inverse multiquadric'
        if length(TYPE) == 1
            C = 1;
        else
            C = TYPE{2};
        end
        
        validateattributes(C,{'double'},{'scalar','nonempty','nonnan',...
            'finite','real','positive'},'computeRBFInterpMatrix',...
            'the value of C when TYPE is ''Inverse Multiquadric''',2);
        A = 1./sqrt((dMat.^2)+C^2);
    case 'multiquadric'
        if length(TYPE) == 1
            C = 1;
        else
            C = TYPE{2};
        end
        
        validateattributes(C,{'double'},{'scalar','nonempty','nonnan',...
            'finite','real','positive'},'computeRBFInterpMatrix',...
            'the value of C when TYPE is ''Inverse Multiquadric''',2);
        A = sqrt((dMat.^2)+C^2);
    case 'linear'
        A = dMat;
    case 'tps'
        A = (dMat.^2).*log(dMat);
    case 'cubic'
        A = dMat.^3;
    otherwise
        error('Unrecognized first entry in TYPE.')
end

end

function inputs = parseComputeRBFInterpMatrixInputs(X,opts)
%PARSECOMPUTERBFINTERPMATRIXINPUTS Parse and verify inputs for the 
%computeRBFInterpMatrix function.

p = inputParser;

% Required inputs
addRequired(p,'X',@(x)validateattributes(X,{'double'},{'2d','nonempty',...
    'nonnan','finite','real'},'computeRBFInterpMatrix','X',1));
X = shiftdim(X); % If X is a vector, force it to be a column vector.

% Optional inputs
addOptional(p,'TYPE',{'Gaussian',20},@(x)validateattributes(x,{'cell'},...
    {'vector','nonempty'},'computeRBFInterpMatrix','TYPE',2));

p.CaseSensitive = false;
p.FunctionName = 'computeRBFInterpMatrix';

parse(p,X,opts{:});

inputs = p.Results;

end
