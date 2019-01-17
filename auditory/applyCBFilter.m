function out = applyCBFilter(in,fS,varargin)
%APPLYCBFILTER Filter input signal using auditory critical-band filters.
%   B = APPLYCBFILTER(A,FS) filters the input signal(s), A, by a series of 
%   Patterson's auditory band filters centered at frequencies between 200 
%   Hz and 16 kHz in intervals of 1 ERB (equivalent rectangular bandwidth).
%   The sampling rate, FS, must also be specified in Hz. 
%   
%   If A is a vector, B is returned as a matrix with the columns containing
%   the input, A, filtered by a Patterson filter at each center frequency.
%
%   If A is a matrix, each signal must be stored as a column, and the 
%   filtering is performed independently on each column. The output, B, is 
%   a cell array, with each element containing a matrix of filtered signals 
%   where each column corresponds to a different center frequency. The 
%   number of elements in the cell array is equal to the number of input
%   signal(s) in A.
%
%   B = APPLYCBFILTER(...,FC) optionally specifies a vector of center 
%   frequencies in Hz.
%
%   B = APPLYCBFILTER(...,FC,N) optionally specifies the length, N, of the
%   auditory band filters to use. N must be specified in seconds. The
%   default length of the filters is equal to the length of the input
%   signal(s) in A.
%
%   B = APPLYCBFILTER(...,FTYPE) optionally specifies the type of auditory
%   band filter to use. FTYPE can take the following inputs:
%       1. 'Patterson' - Patterson's auditory filters (default)
%       2. 'Gammatone' - Gammatone filters (requires LTFAT toolbox)
%
%   B = APPLYCBFILTER(...,'sequential') sequentially filters the input
%   signal through each critical band. This produces an output signal, B,
%   that has the same dimensions as the input signal(s), A.
%
%   See also GETGAMMATONEFILTERS, GETPATTERSONFILTERS,
%   APPLYFILTERSINSERIES.

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

narginchk(2,5);

% Parse and verify inputs
[inputs,extra] = parseApplyCBFilterInputs(in,fS,varargin);

% Extract parsed inputs
in = inputs.in;
fS = inputs.fS;
fC = inputs.fC;
n = inputs.n;
fType = inputs.fType;
METHOD = inputs.METHOD;

if isempty(fC)
    fC = extra{1,1};
end

filtLen = round(n*fS);
switch lower(fType)
    case 'patterson'
        filtMat = getPattersonFilters(fC,fS,filtLen);
    case 'gammatone'
        filtMat = getGammatoneFilters(fC,fS,filtLen);
    otherwise
        error('Unrecognized input for FTYPE.')
end

[sigLen,numCh] = size(in);
numFCs = size(filtMat,2);
switch lower(METHOD)
    case 'default'
        out = cell(numCh,1);
        for ii = 1:numCh
            out{ii,1} = zeros(sigLen,numFCs);
            for jj = 1:numFCs
                out{ii,1}(:,jj) = filter(filtMat(:,jj),1,in(:,ii));
            end
        end
        
        if numCh == 1 % For a single input signal, return a matrix.
            out = cell2mat(out);
        end
    case 'sequential'
        out = applyFiltersInSeries(filtMat,in);
    otherwise
        error('Unrecognized input.')
end

end

function [inputs,extra] = parseApplyCBFilterInputs(in,fS,opts)
%PARSEAPPLYCBFILTERINPUTS Parse and verify inputs for the applyCBFilter 
%function.

p = inputParser;

% Required inputs
addRequired(p,'in',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite','real'},'applyCBFilter','A',1));
addRequired(p,'fS',@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','nonnan','finite','real','positive'},'applyCBFilter',...
    'FS',2));

% If 'in' is a row vector, force it to be a column vector.
in = shiftdim(in);

% Check if fS is specified in Hz.
if fS < 1000
    warning(['Low value of FS detected. Note that FS must be specified',...
        ' in Hz.'])
end

% Optional inputs
fCVec = getERBFreqVec(200,16000);
addOptional(p,'fC',fCVec,@(x)validateattributes(x,{'double'},{'vector',...
    'nonnan','finite','real','nonnegative'},'applyCBFilter','FC',3));
irLen = size(in,1);
addOptional(p,'n',irLen/fS,@(x)validateattributes(x,{'double'},{...
    'scalar','nonnan','finite','real','positive'},'applyCBFilter','N',4));
addOptional(p,'fType','Patterson',@(x)validateattributes(x,{'char'},{...
    'scalartext','nonempty'},'applyCBFilter','FTYPE',5));
addParameter(p,'METHOD','default',@(x)validateattributes(x,{'char'},{...
    'scalartext','nonempty'},'applyCBFilter'));

% Return additional variables that may be used in main function
extra{1,1} = fCVec;

p.CaseSensitive = false;
p.FunctionName = 'applyCBFilter';

parse(p,in,fS,opts{:});

inputs = p.Results;

end
