function [retFlag, retMsg, retDat] = sendOSCmessage(udps, sendMsg, sendDat, udpr, retMsg, retDat, waitTime)
%SENDOSCMESSAGE Send an OSC message.
%   FLAG = SENDOSCMESSAGE(UDPS, SMSG, SDAT) sends OSC message and data
%   (SMSG, SDAT) over UDPSender UDPS. Returns FLAG as false since to
%   indicate that no response is received (see below).
%
%   FLAG = SENDOSCMESSAGE(UDPS, SMSG, SDAT, UDPR, RMSG) waits for up to 0.5
%   seconds for a response message (RMSG) to arrive over UDPReceiver UDPR.
%   FLAG is true only if the desired response message is received within
%   the wait time.
%
%   FLAG = SENDOSCMESSAGE(UDPS, SMSG, SDAT, UDPR, RMSG, RDAT) waits for
%   the specified response message and data (RMSG, RDAT) to arrive. FLAG is
%   true only if the desired message and data are received together within
%   the wait time.
%
%   FLAG = SENDOSCMESSAGE(UDPS, SMSG, SDAT, UDPR, RMSG, RDAT, T) waits for
%   T seconds for the desired response instead of the usual 0.5 seconds.
%
%   [FLAG, RMSG, RDAT] = SENDOSCMESSAGE(UDPS, SMSG, SDAT, UDPR, RMSG, ...)
%   returns also the received OSC message and data (RMSG, RDAT).
%
%   See also WRITEOSCMESSAGE, READOSCMESSAGE.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2017 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   =======================================================================

if nargin < 7
    waitTime = 0.5;
end

if nargin < 6
    retDat = [];
end

if nargin < 5
    retMsg = [];
end

sendBuff = writeOSCmessage(sendMsg, sendDat);
step(udps,sendBuff);

retFlag = false;
if ~isempty(retMsg)
    % Wait for return message
    tic;
    while toc < waitTime
        retBuff = step(udpr);
        if ~isempty(retBuff)
            if ~isempty(retDat)
                % Check return message and data
                [msg, dat] = readOSCmessage(retBuff);
                if strcmpi(msg,retMsg) && isequal(retDat, dat)
                    retFlag = true;
                    break;
                end
            else
                % Check return message only
                [msg, retDat] = readOSCmessage(retBuff);
                if strcmpi(msg,retMsg)
                    retFlag = true;
                    break;
                end
            end
        end
    end
    
    if ~retFlag
        warning('Return message ''%s'' was not received!', retMsg);
    end
end

end