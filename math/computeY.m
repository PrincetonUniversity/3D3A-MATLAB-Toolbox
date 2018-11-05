function Y = computeY(n,dirs,varargin)
%COMPUTEY Compute spherical harmonics.
%   Y = COMPUTEY(n,dirs) returns real-valued spherical harmonics of degree
%   n for the directions specified in dirs. n may be a vector but each of 
%   its elements must be non-negative. Y is a cell array with length equal 
%   to the length of n. Y{1,ii} contains the spherical harmonics of degree 
%   n(ii), where ii is an indexing variable. dirs must be a p-by-3 matrix 
%   of directions specified in SOFA cartesian coordinates. Each Y{1,ii} 
%   contains a matrix of dimensions p-by-(2*n(ii)+1). The returned 
%   spherical harmonics are orthonormal. The same command may also be 
%   specified as: 
%       Y = COMPUTEY(n,dirs,'TYPE','real').
%
%   Y = COMPUTEY(...,'TYPE','complex') returns complex-valued spherical
%   harmonics instead.
%
%   Y = COMPUTEY(...,'CSPHASE',0) ignores the Condon-Shortley phase term 
%   (default). 
%
%   Y = COMPUTEY(...,'CSPHASE',1) includes the Condon-Shortley phase term.
%
%   See also COMPUTEYMAT.

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
inputs = parseCOMPUTEYInputs(n,dirs,varargin);

% Extract parsed inputs
n = shiftdim(inputs.n); % Force n to be a column vector
dirs = inputs.dirs;
TYPE = inputs.TYPE;
CSPHASE = inputs.CSPHASE;

% Main computation begins

numDegrees = length(n);
numDirs = size(dirs,1);
[az,el,~] = cart2sph(dirs(:,1),dirs(:,2),dirs(:,3));

Y = cell(1,numDegrees);
if strcmpi(TYPE,'complex')
    for ii = 1:numDegrees
        numOrders = (2*n(ii))+1;
        m = (linspace(-n(ii),n(ii),numOrders)).';
        normTerm = sqrt((numOrders/(4*pi))*(factorial(n(ii)-abs(m)))./...
            factorial(n(ii)+abs(m)));
        signTerm = (-1).^abs(m); % Cancel CS phase term in legendre
        if CSPHASE == 1
            signTerm(ceil(numOrders/2):end) = 1; % Don't cancel CS phase
        end
        Pnm = zeros(numOrders,numDirs);
        Pnm(ceil(numOrders/2):end,:) = legendre(n(ii),sin(el));
        % The legendre function computes associated Legendre functions, 
        % nPm, for positive values of m for a given n. Since, for negative 
        % m, we have nP(-m) = signTerm*normTerm*nPm, we need only flip the 
        % above calculations to populate the first half of Pnm below.
        if numOrders > 1
            Pnm(1:floor(numOrders/2),:) = flipud(Pnm(...
                ceil(numOrders/2)+1:end,:));
        end
        phaseTerm = exp(1i*m*(az)');
        Y{1,ii} = (repmat(signTerm,1,numDirs).*repmat(normTerm,1,...
            numDirs).*Pnm.*phaseTerm).';
    end
elseif strcmpi(TYPE,'real')
    for ii = 1:numDegrees
        numOrders = (2*n(ii))+1;
        m = (linspace(-n(ii),n(ii),numOrders)).';
        normTerm = sqrt((numOrders*(2-(~m))/(4*pi)).*...
            (factorial(n(ii)-abs(m)))./factorial(n(ii)+abs(m)));
        signTerm = (-1).^abs(m); % Cancel CS phase term in legendre
        if CSPHASE == 1 
            signTerm(ceil(numOrders/2):end) = 1; % Don't cancel CS phase
        end
        Pnm = zeros(numOrders,numDirs);
        Tnm = zeros(numOrders,numDirs);
        Pnm(ceil(numOrders/2):end,:) = legendre(n(ii),sin(el));
        Tnm(ceil(numOrders/2):end,:) = cos(m(ceil(numOrders/2):end)...
            *((az).'));
        % See note on legendre function calculation for the TYPE =
        % 'complex' case above for an explanation of the following two
        % calculations.
        if numOrders > 1
            Pnm(1:floor(numOrders/2),:) = flipud(Pnm(...
                ceil(numOrders/2)+1:end,:));
            Tnm(1:floor(numOrders/2),:) = ...
                sin(abs(m(1:floor(numOrders/2)))*((az).'));
        end
        Y{1,ii} = (repmat(signTerm,1,numDirs).*repmat(normTerm,1,...
            numDirs).*Pnm.*Tnm).';
    end
else
    error('Unrecognized input. TYPE must be ''%s'' or ''%s.''',...
        'real','complex')
end

% Main computation ends

end

function inputs = parseCOMPUTEYInputs(n,dirs,opts)
%PARSECOMPUTEYINPUTS Parse and verify inputs for the computeY function.

p = inputParser;

% Required inputs
addRequired(p,'n',@(x)validateattributes(x,{'double'},{'vector',...
    'nonempty','nonnan','finite','nonnegative'},'computeY','n',1));
addRequired(p,'dirs',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite','size',[NaN,3]},'computeY','dirs',2));

% Optional inputs

addParameter(p,'TYPE','real',@(x)validateattributes(x,{'char'},...
    {'nonempty'},'computeY','TYPE'));
addParameter(p,'CSPHASE',0,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','integer','nonnegative','<=',1},'computeY',...
    'CSPHASE'));

p.CaseSensitive = false;
p.FunctionName = 'computeY';

parse(p,n,dirs,opts{:});

inputs = p.Results;

end
