function [outputData,YCoeffs,YMat,degVal] = getYRep(inputData,Ydata,...
    varargin)
%GETYREP Compute a spherical harmonic representation of input data.
%   [outputData,YCoeffs,YMat,degVal] = GETYREP(inputData,{'Ys',YMat}) 
%   returns the spherical harmonic representation of the data in inputData 
%   using the matrix of spherical harmonics in YMat. If inputData has 
%   dimensions N-by-p, YMat must have dimensions p-by-M. See COMPUTEYMAT 
%   for a function to compute YMat. Also returned are the spherical 
%   harmonic coefficients, YCoeffs. outputData will have the same 
%   dimensions as inputData, while YCoeffs will have dimensions M-by-N. The 
%   spherical harmonic degree used to represent the input data is returned 
%   as degVal (useful when the 'checkMaxD' option is set to 1). The 
%   following optional input may also be provided:
%       GETYREP(...,'checkMaxD',1) - limit the max. spherical harmonic
%       degree used to ensure that the calculation of YCoeffs is an
%       overdetermined problem. This check is not performed by default.
%
%   ___ = GETYREP(inputData,{'YDirs',inputDirs,'D',maxD}) returns the maxD-
%   degree spherical harmonic representation of the data in inputData whose 
%   values are specified for the spatial directions given in inputDirs. If 
%   inputData has dimensions N-by-p, inputDirs must be specified in SOFA 
%   cartesian coordinates as a p-by-3 matrix, where p is the number of 
%   directions. maxD must be a non-negative scalar. In addition to 
%   returning outputData and YCoeffs, the matrix of spherical harmonics, 
%   YMat, computed internally using inputDirs and maxD, are also returned. 
%   YMat will have dimensions p-by-(maxD+1)^2. The following optional 
%   inputs may also be provided:
%       1. GETYREP(...,'TYPE','real') - use real-valued spherical
%       harmonics (default).
%       2. GETYREP(...,'TYPE','complex') - use complex-valued spherical
%       harmonics.
%       3. GETYREP(...,'CSPHASE',0) - Ignore Condon-Shortley phase term
%       (default).
%       4. GETYREP(...,'CSPHASE',1) - Include Condon-Shortley phase term.
%       5. GETYREP(...,'checkMaxD',1) - limit the max. spherical harmonic
%       degree used to ensure that the calculation of YCoeffs is an
%       overdetermined problem. This check is not performed by default.
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

narginchk(2,8);

% Parse and verify inputs
inputs = parseGETYREPInputs(inputData,Ydata,varargin);

% Extract parsed inputs
inputData = inputs.inputData;
Ydata = inputs.Ydata;
typeVal = inputs.TYPE;
csPhaseFlag = inputs.CSPHASE;
checkMaxDFlag = inputs.checkMaxD;

nCols = size(inputData,2);
switch length(Ydata)
    case 2
        validateattributes(Ydata{1},{'char'},{'scalartext'},'getYRep',...
            'Ydata{1}')
        if strcmpi(Ydata{1},'Ys')
            validateattributes(Ydata{2},{'numeric'},{'2d','nonempty',...
                'nonnan','finite','size',[nCols,NaN]},'getYRep','Ydata{2}')
            YMat = Ydata{2};
            if ~isempty(varargin)
                warning(['Ignoring all but first 2 inputs since Ydata',...
                    ' is specified in the format: {''Ys'',YMat}.'])
            end
        else
            error(['Unrecognized input. First element in Ydata ',...
                'must be ''%s'' when Ydata has length 2.'],'Ys')
        end
        
        inferredD = sqrt(size(YMat,2))-1; % Infer the max. spherical 
                                          % harmonic degree from YMat.
        if floor(inferredD) == inferredD % Check if inferredD is an integer
            
            % When checkMaxDFlag is true, the spherical harmonic degree is
            % limited to ensure that solving for YCoeffs later is an 
            % overdetermined problem.
            if checkMaxDFlag
                maxD = floor(sqrt(nCols/2))-1;
                inferredD = min([inferredD,maxD]);
            end
            degVal = inferredD;
        else
            error(['Number of columns in YMat must be a perfect ',...
                'square.']);
        end
    case 4
        validateattributes(Ydata{1},{'char'},{'scalartext'},'getYRep',...
            'Ydata{1}')
        if strcmpi(Ydata{1},'YDirs')
            validateattributes(Ydata{2},{'numeric'},{'2d','nonempty',...
                'nonnan','finite','size',[nCols,3]},'getYRep','Ydata{2}')
            inputDirs = Ydata{2};
        else
            error(['Unrecognized input. First element in Ydata ',...
                'must be ''%s'' when Ydata has length 4.'],'YDirs')
        end
        validateattributes(Ydata{3},{'char'},{'scalartext'},'getYRep',...
            'Ydata{3}')
        if strcmpi(Ydata{3},'D')
            validateattributes(Ydata{4},{'numeric'},{'scalar',...
                'nonempty','nonnan','finite','nonnegative'},'getYRep',...
                'Ydata{4}')
            maxD = Ydata{4};
        else
            error(['Unrecognized input. Third element in Ydata ',...
                'must be ''%s'' when Ydata has length 4.'],'D')
        end
        if strcmpi(typeVal,'real') || strcmpi(typeVal,'complex')
            
            % See checkMaxDFlag comment under case 2.
            if checkMaxDFlag
                maxValidD = floor(sqrt(nCols/2))-1;
                validD = min([maxD,maxValidD]);
            else
                validD = maxD;
            end
            degVal = validD;
            
            % Compute matrix of spherical harmonics up to degree maxD
            % sampled at directions specified in inputDirs.
            YMat = computeYMat(degVal,inputDirs,'TYPE',typeVal,...
                'CSPHASE',csPhaseFlag);
        else
            error('Unrecognized input. TYPE must be ''%s'' or ''%s.''',...
                'real','complex')
        end
    otherwise
        error('Invalid dimensions for Ydata.')
end

% Main computation begins

validDs = 1:((degVal+1)^2);
YMat = YMat(:,validDs);
invYMat = pinv(YMat);
YCoeffs = invYMat*inputData.';
outputData = (YMat*YCoeffs).';

% Main computation ends

end

function inputs = parseGETYREPInputs(inputData,Ydata,opts)
%PARSEGETYREPINPUTS Parse and verify inputs for the getYRep function.

p = inputParser;

% Required inputs
addRequired(p,'inputData',@(x)validateattributes(x,{'double'},{'2d',...
    'nonempty','nonnan','finite'},'getYRep','inputData',1));
addRequired(p,'Ydata',@(x)validateattributes(x,{'cell'},{'nonempty'},...
    'getYRep','Ydata',2));

% Optional inputs
if length(Ydata) == 4
    addParameter(p,'TYPE','real',@(x)validateattributes(x,{'char'},...
    {'nonempty'},'getYRep','TYPE'));
    addParameter(p,'CSPHASE',0,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','integer','nonnegative','<=',1},'getYRep',...
    'CSPHASE'));
    addParameter(p,'checkMaxD',0,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','integer','nonnegative','<=',1},'getYRep',...
    'checkMaxD'));
else
    addParameter(p,'TYPE','real'); % Not used, so no validation performed.
    addParameter(p,'CSPHASE',0); % Not used, so no validation performed.
    addParameter(p,'checkMaxD',0); % Not used, so no validation performed.
end

p.CaseSensitive = false;
p.FunctionName = 'getYRep';

parse(p,inputData,Ydata,opts{:});

inputs = p.Results;

end
