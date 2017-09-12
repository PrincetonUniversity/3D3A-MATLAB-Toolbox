function buff = writeOSCmessage(msg, dat)
%writeOSCmessage Encode an OSC message.
%   B = writeOSCmessage(MSG, DAT) encodes the OSC message (MSG) and data
%   (DAT) into a character buffer (B) for sending across UDP.
%
%   See also readOSCmessage.

oscAddressPattern = uint8(makeOSCstring(msg));

if ~iscell(dat)
    dat = {dat};
end
numArgs = length(dat);

oscTypeTags = uint8(zeros(1,ceil((1+numArgs)/4)*4));
oscArguments = cell(1,numArgs);

oscTypeTags(1) = uint8(',');
for ii = 1:numArgs
    if isinteger(dat{ii})
        % Integer data
        oscTypeTags(1+ii) = uint8('i');
        oscArguments{ii} = flip(typecast(dat{ii},'uint8'));
    elseif isfloat(dat{ii})
        % Float data
        oscTypeTags(1+ii) = uint8('f');
        oscArguments{ii} = flip(typecast(single(dat{ii}),'uint8'));
    elseif ischar(dat{ii})
        % String data
        oscTypeTags(1+ii) = uint8('s');
        oscArguments{ii} = uint8(makeOSCstring(dat{ii}));
    else
        warning('Unknown OSC message data type.');
    end
end

buff = [oscAddressPattern oscTypeTags cell2mat(oscArguments)].';

end

function s = makeOSCstring(s)
s = [s 0 0 0 0];
s = s(1:end-mod(length(s),4));
end