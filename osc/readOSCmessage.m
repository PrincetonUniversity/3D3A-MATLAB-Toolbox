function [addr, args] = readOSCmessage(buff)
%readOSCmessage Parse an OSC message.
%   [MSG, DAT] = readOSCmessage(B) parses the OSC character buffer (B) into
%   the OSC message (MSG) and the corresponding data (DAT).
%
%   See also writeOSCmessage.

msgLen = length(buff);
charmsg = char(buff.');
commaIndex = find(charmsg==',',1,'first');
oscAddressPatternEnd = find(buff(1:commaIndex-1),1,'last');
oscAddressPattern = charmsg(1:oscAddressPatternEnd);

oscTypeTagStringEnd = commaIndex + 3;
while(buff(oscTypeTagStringEnd)~=0)
    oscTypeTagStringEnd = oscTypeTagStringEnd + 4;
end
numArgs = find(buff(commaIndex:oscTypeTagStringEnd),1,'last') - 1;
oscTypeTagString = charmsg(commaIndex:commaIndex + numArgs);

oscArguments = cell(1,numArgs);
argStart = oscTypeTagStringEnd + 1;
for ii=1:numArgs
    argEnd = argStart + 3;
    if(argEnd > msgLen)
        disp('Warning: Attempted to exceed OSC message length.')
        break
    end
    switch oscTypeTagString(ii + 1)
        case 'i'
            % Integer data
            oscArguments{ii} = typecast(flip(buff(argStart:argEnd)),'int32');
        case 'f'
            % Float data
            oscArguments{ii} = typecast(flip(buff(argStart:argEnd)),'single');
        case 's'
            % String data
            while(buff(argEnd)~=0)
                argEnd = argEnd + 4;
            end
            strEnd = argStart + find(buff(argStart:argEnd),1,'last') - 1;
            oscArguments{ii} = charmsg(argStart:strEnd);
        otherwise
            warning('Unknown OSC message data type.');
            oscArguments{ii} = 0;
    end
    argStart = argEnd + 1;
end

addr = oscAddressPattern;
args = oscArguments;

end