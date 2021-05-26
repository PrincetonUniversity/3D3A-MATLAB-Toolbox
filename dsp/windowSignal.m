function [hWin,varargout] = windowSignal(h,wLen,varargin)
%WINDOWSIGNAL Window a signal.
%   B = WINDOWSIGNAL(A,L) applies a rectangular window of length L samples 
%   to the signal(s) in A. The window is applied starting at sample 1 for 
%   all signals in A.
%       If A is a row or column vector, B will also be a row or column 
%       vector, respectively.
%       If A is a matrix, the individual signals must be stored as columns. 
%       B will then contain the windowed signals stored as column vectors.
%
%   The length of the signals in B are equal to L. If L exceeds the length 
%   of the signals in A, the signals in A are first windowed and then zero-
%   padded on the right.
%
%   B = WINDOWSIGNAL(...,'wType',TYPE) optionally specifies the type of
%   window to apply. The following may be specified for TYPE (which must be
%   a cell array):
%       1. {'rect'} - rectangular window (default)
%       2. {'hann'} - hann window
%       3. {'hamming'} - hamming window
%       4. {'tukey',R} - tukey window with the first and last 100*R/2 
%       percent of samples equal to parts of a cosine. If R is not 
%       specified, a default of 0.5 is assumed. R can take any value in the 
%       range 0 to 1. For more information, see TUKEYWIN.
%       5. {'rc',[R1,R2]} - raised-cosine window with the first 100*R1
%       and last 100*R2 percent of samples equal to parts of a cosine. If
%       R1 and R2 are not specified, defaults of R1 = 0.25 and R2 = 0.25
%       are assumed, which corresponds to a Tukey window with R = 0.5. R1
%       and R2 can take values in the range 0 to 1 such that R1+R2 <= 1. 
%       For more information, see RAISEDCOSINEWIN.
%
%   B = WINDOWSIGNAL(...,'start',S) optionally specifies the sample number
%   at which the window is to be applied to the signal(s) in A. S must be a 
%   positive integer or a vector of positive integers with length equal to
%   the number of columns in A.
%
%   B = WINDOWSIGNAL(...,'offset',O) optionally specifies additional sample
%   amounts by which the window is offset prior to being applied to the 
%   signal(s) in A. O must be a vector with length equal to the number of 
%   columns in A. The elements of O must be non-negative integers. The
%   smallest value of the specified O must be 0. If not, the specified O
%   will be rescaled in this way.
%
%   [B,W] = WINDOWSIGNAL additionally returns W, the matrix of 
%   appropriately positioned windows used to window each of the signals in 
%   A. W will have the same dimensions as B.
%
%   EXAMPLE: 
%
%       Suppose a pair of impulse responses, IR1 and IR2, have length N 
%       samples each. Their onsets occur at sample numbers M1 and M2, 
%       respectively, where M2-M1 = P > 0. Therefore, the "tail" (i.e., 
%       which we use to refer to the portion of the IR after the onset) of 
%       IR2 is shorter than that of IR1 by P samples, and each IR has an 
%       effective length of N-M1 (for IR1) and N-M2 (for IR2). To window 
%       both IRs to have a new length of Q, it may be desirable to apply 
%       the window such that the following conditions are satisfied:
%           (i) The window begins R < min([M1,M2]) samples before the main 
%           impulse for BOTH IR1 and IR2.
%           (ii) The relative delay of P samples between IR1 and IR2 must
%           be preserved.
%       
%       To window the two IRs to meet both conditions above, specify the
%       call to this function as follows (assuming IR1 and IR2 are each 
%       column vectors, and the window is a raised-cosine window with 
%       default window parameters):
%           [B,W] = windowSignal([IR1,IR2],Q,'wType',{'rc'},'start',...
%               M1-R+1,'offset',[0,P]);
%
%       To window the two IRs while meeting only condition (i) above,
%       specify the call to this function as follows:
%           [B,W] = windowSignal([IR1,IR2],Q,'wType',{'rc'},'start',...
%               [M1-R+1,M2-R+1]);
%       In this case, the pre-onset delay and the lengths of the "tails"
%       for IR1 and IR2, after windowing, will be the same but the relative
%       delay, P = M2-M1, between the IRs will be lost.
%
%       To window the two IRs while meeting only condition (ii) above, 
%       specify the call to this function as follows:
%           [B,W] = windowSignal([IR1,IR2],Q,'wType',{'rc'},'start',...
%               M1-R+1);
%       In this case, the pre-onset delay and the lengths of the "tails" 
%       for IR1 and IR2, after windowing, will be different.
%
%   Needs: Signal Processing Toolbox.
%
%   See also HANN, HAMMING, TUKEYWIN, RAISEDCOSINEWIN.

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
%   Copyright (c) 2021 Princeton University
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

% Validate required inputs
validateattributes(h,{'numeric'},{'2d','nonempty','nonnan','finite'},...
    'windowSignal','A',1)
validateattributes(wLen,{'numeric'},{'scalar','finite','positive',...
    'integer'},'windowSignal','L',2)

if isrow(h)
    h = h.';
    transposeFlag = true;
else
    transposeFlag = false;
end

[hLen,numCols] = size(h);

% Validate optional inputs
indx = find(strcmpi(varargin,'wType'),1);
if ~isempty(indx)
    wType = varargin{indx+1};
    validateattributes(wType,{'cell'},{'nonempty','vector'},...
        'windowSignal','TYPE for option: wType')
else
    wType = {'rect'};
end

indx = find(strcmpi(varargin,'start'),1);
if ~isempty(indx)
    wS = varargin{indx+1};
    validateattributes(wS,{'numeric'},{'vector','finite','positive',...
        'integer'},'windowSignal','S for option: start')
else
    wS = ones(1,numCols);
end

indx = find(strcmpi(varargin,'offset'),1);
if ~isempty(indx)
    oVec = varargin{indx+1};
    validateattributes(oVec,{'numeric'},{'vector','finite',...
        'nonnegative','integer','numel',numCols},'windowSignal',...
        'O for option: offset')
else
    oVec = zeros(1,numCols);
end

if isscalar(wS)
    wS = wS*ones(1,numCols);
end

% Rescale oVec so that the smallest value is 0.
oVec = oVec-min(oVec);

% Compute longest admissable window length 
tS = wS+oVec; % Desired starting sample for each windowed IR
maxStartSample = max(tS);
maxWinLen = hLen-maxStartSample+1;
appWinLen = min([wLen,maxWinLen]); % Longest admissable window length

% Create window vector
validateattributes(wType{1,1},{'char'},{'scalartext','nonempty'},...
    'windowSignal','TYPE{1} for option: wType')
switch lower(wType{1,1})
    case 'rect'
        wVec = rectwin(appWinLen); % From Signal Processing Toolbox
    case 'hann'
        wVec = hann(appWinLen); % From Signal Processing Toolbox
    case 'hamming'
        wVec = hamming(appWinLen); % From Signal Processing Toolbox
    case 'tukey'
        if length(wType) > 1
            validateattributes(wType{1,2},{'numeric'},{'scalar','real',...
                'nonnegative','<=',1},'windowSignal',['TYPE{2} when ',...
                'TYPE{1} = ''tukey'' for option: wType']);
            R = wType{1,2};
        else
            R = 0.5;
        end
        wVec = tukeywin(appWinLen,R); % From Signal Processing Toolbox
    case 'rc'
        if length(wType) > 1
            validateattributes(wType{1,2},{'numeric'},{'vector','real',...
                'nonnegative','<=',1,'numel',2},'windowSignal',...
                'TYPE{2} when TYPE{1} = ''rc'' for option: wType');
            R1 = wType{1,2}(1);
            R2 = wType{1,2}(2);
        else
            R1 = 0.25;
            R2 = 0.25;
        end
        wVec = raisedcosinewin(appWinLen,[R1,R2]);
    otherwise
        error('Invalid TYPE specification for option: wType.')
end

% Apply window
hWin = zeros(wLen,numCols);
winMat = zeros(wLen,numCols);
for ii = 1:numCols
   winStartIndx = oVec(ii)+1;
   hWin(winStartIndx:winStartIndx+appWinLen-1,ii) = h(tS(ii):tS(ii)+...
       appWinLen-1,ii).*wVec;
   winMat(winStartIndx:winStartIndx+appWinLen-1,ii) = wVec;
end

if transposeFlag
    hWin = hWin.';
    winMat = winMat.';
end

if nargout > 1
    varargout{1} = winMat;
end

end
