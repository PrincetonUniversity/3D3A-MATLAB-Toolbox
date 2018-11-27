function [hWin,varargout] = windowSignal(h,wLen,varargin)
%WINDOWSIGNAL Window a signal.
%   hWin = WINDOWSIGNAL(h,wLen) applies a rectangular window of length wLen 
%   (specified in samples) to the signal(s) in h. The window is applied  
%   starting at sample 1 for all signals in h. If h is a matrix, the 
%   individual signals must be stored as columns. hWin contains the 
%   windowed signals stored as column vectors. The length of the signals in 
%   hWin are equal to wLen.
%
%   ___ = WINDOWSIGNAL(...,'wType',TYPE) optionally specifies the type of
%   window to apply. The following may be specified for TYPE (which must be
%   a cell array):
%       1. {'rect'} - rectangular window (default)
%       2. {'hann'} - hanning window
%       3. {'hamm'} - hamming window
%       4. {'tukey',R} - tukey window with the first 100*R/2 percent of
%       samples equal to parts of a cosine. If R is not specified, a
%       default of 0.5 is assumed. R can take values in the range 0 to 1. 
%       For more information, see TUKEYWIN.
%       5. {'rc',[R1,R2]} - raised-cosine window with the first 100*R1
%       and last 100*R2 percent of samples equal to parts of a cosine. If
%       R1 and R2 are not specified, defaults of R1 = 0.25 and R2 = 0.25
%       are assumed, which corresponds to a Tukey window with R = 0.5. R1
%       and R2 can take values in the range 0 to 0.5. For more information,
%       see RAISEDCOSINEWIN.
%
%   ___ = WINDOWSIGNAL(...,'start',wS) optionally specifies the sample 
%   number at which the window is to be applied to the signal(s) in h. wS
%   must be a positive integer.
%
%   ___ = WINDOWSIGNAL(...,'offset',oVec) optionally specifies additional
%   sample amounts by which the window is offset prior to being applied to 
%   the signal(s) in h. oVec must be a vector with length equal to the 
%   number of columns in h. The elements of oVec must be non-negative
%   integers.
%
%   [___,L] = WINDOWSIGNAL optionally returns the actual length of the
%   window used to window the signals in h. L will equal wLen if oVec is
%   not specified or contains all zeros. Otherwise, L will be smaller than
%   wLen by an amount equal to max(oVec).
%
%   EXAMPLE: 
%
%       Suppose a pair of impulse responses, IR1 and IR2, have length N 
%       samples each. The main impulse in IR1 and IR2 occur at sample 
%       numbers M1 and M2, respectively, where M2-M1 = P > 0. Therefore, 
%       the "tail" of IR2 is shorter than that of IR1 by P samples, and 
%       each IR has an effective length of N-M1 (for IR1) and N-M2 (for 
%       IR2). To window both IRs to have a new length of Q, where Q < 
%       (N-M2), it may be desirable to apply the window such that the 
%       following conditions are satisfied:
%           (i) The window begins R < min([M1,M2]) samples before the main 
%           impulse for both IR1 and IR2. Note: a naive application of such
%           a window will eliminate any relative delay between the two IRs.
%           (ii) The relative delay of P samples between IR1 and IR2 must
%           be preserved.
%       
%       To window the two IRs to meet both conditions above, specify the
%       call to the WINDOWSIGNAL function as follows (assuming IR1 and IR2
%       are each column vectors, and the window is a raised-cosine window 
%       with default window parameters):
%           [hWin,L] = windowSignal([IR1,IR2],Q,'wType',{'rc'},'start',...
%               M1-R+1,'offset',[0,P]);
%       In this case, the effective window length, L = Q-P.
%
%       To window the two IRs while meeting only condition (ii) above and
%       requiring L = Q, specify the call to the WINDOWSIGNAL function as
%       follows:
%           hWin = windowSignal([IR1,IR2],Q,'wType',{'rc'},'start',M1-R+1);
%
%   See also TUKEYWIN, RAISEDCOSINEWIN.

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
inputs = parseWindowSignalInputs(h,wLen,varargin);

% Extract parsed inputs
h = inputs.h;
wLen = inputs.wLen;
wType = inputs.wType;
wS = inputs.start;
oVec = inputs.offset;

h = shiftdim(h); % If h is a vector, force to a column.
[hLen,numCols] = size(h);
oVecMax = max(oVec);
if wLen > (hLen-wS-oVecMax+1)
    error(['The specified wLen is %d samples. The maximum allowable ',...
        'window length is length(h)-wS-max(oVec)+1 = %d samples for ',...
        'the provided inputs.'],wLen,hLen-wS-oVecMax+1)
end
validWLen = wLen-oVecMax;

% Create window vector
validateattributes(wType{1},{'char'},{'scalartext','nonempty'},...
    'windowSignal','TYPE{1} for option: wType')
switch lower(wType{1})
    case 'rect'
        wVec = ones(validWLen,1);
    case 'hann'
        wVec = hann(validWLen);
    case 'hamm'
        wVec = hamming(validWLen);
    case 'tukey'
        if length(wType) > 1
            validateattributes(wType{2},{'double'},{'scalar','nonempty',...
                'real','nonnegative','<=',1},'windowSignal',['TYPE{2} ',...
                'when TYPE{1} = ''tukey'' for option: wType']);
            R = wType{2};
        else
            R = 0.5;
        end
        wVec = tukeywin(validWLen,R);
    case 'rc'
        if length(wType) > 1
            validateattributes(wType{2},{'double'},{'vector','nonempty',...
                'real','nonnegative','<=',0.5,'numel',2},'windowSignal',...
                'TYPE{2} ','when TYPE{1} = ''rc'' for option: wType');
            R1 = wType{2}(1);
            R2 = wType{2}(2);
        else
            R1 = 0.25;
            R2 = 0.25;
        end
        wVec = raisedcosinewin(validWLen,[R1,R2]);
    otherwise
        error('Invalid TYPE specification for option: wType.')
end

hWin = zeros(wLen,numCols);
baselineWin = circshift([wVec;zeros(hLen-validWLen,1)],(wS-1));
for ii = 1:numCols
    currentWin = circshift(baselineWin,oVec(ii));
    hWinFull = h(:,ii).*currentWin;
    hWin(:,ii) = hWinFull(wS:(wS+wLen-1));
end

if nargout == 2
    varargout{1} = validWLen;
end

end

function inputs = parseWindowSignalInputs(h,wLen,opts)
%PARSEWINDOWSIGNALINPUTS Parse and verify inputs to the windowSignal
%function.

p = inputParser;

% Required inputs
addRequired(p,'h',@(x)validateattributes(x,{'double'},{'2d','nonempty',...
    'nonnan','finite'},'windowSignal','h',1));
addRequired(p,'wLen',@(x)validateattributes(x,{'double'},{'scalar',...
    'nonempty','nonnan','finite','positive','integer'},'windowSignal',...
    'wLen',2));
numCols = size(shiftdim(h),2);

% Optional inputs
addParameter(p,'wType',{'rect'},@(x)validateattributes(x,{'cell'},...
    {'nonempty','vector'},'windowSignal','TYPE for option: wType'));
addParameter(p,'start',1,@(x)validateattributes(x,{'double'},...
    {'scalar','nonempty','nonnan','finite','positive','integer'},...
    'windowSignal','wS for option: start'));
addParameter(p,'offset',zeros(1,numCols),@(x)validateattributes(x,...
    {'double'},{'vector','nonempty','nonnan','finite','nonnegative',...
    'integer','numel',numCols},'windowSignal','oVec for option: offset'));

p.CaseSensitive = false;
p.FunctionName = 'windowSignal';

parse(p,h,wLen,opts{:});

inputs = p.Results;

end
