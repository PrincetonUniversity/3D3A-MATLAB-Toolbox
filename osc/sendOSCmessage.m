function [retFlag, retMsg, retDat] = sendOSCmessage(udps, sendMsg, sendDat, udpr, retMsg, retDat, waitTime)
%sendOSCmessage Closed-loop OSC messaging.
%   [FLAG, RMSG, RDAT] = sendOSCmessage(UDPS, SMSG, SDAT, UDPR, RMSG, RDAT, T)
%   sends OSC message and data (SMSG, SDAT) over UDPSender UDPS and
%   awaits a return message (RMSG, RDAT) over UDPReceiver UDPR. The maximum
%   wait time for a response is T, whose default value is 0.5 seconds.
%   Returns FLAG as true if the desired message and data are received,
%   otherwise FLAG is false.
%
%   [FLAG, RMSG, RDAT] = sendOSCmessage(UDPS, SMSG, SDAT, UDPR, RMSG, [], T)
%   if the return data is not specified, then any response with the correct
%   return message will be accepted.
%
%   [FLAG, RMSG, RDAT] = sendOSCmessage(UDPS, SMSG, SDAT, UDPR, [], [], T)
%   if neither the return message or data are specified, then this becomes
%   an open-loop messaging function (no response is needed).
%
%   See also writeOSCmessage, readOSCmessage.

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