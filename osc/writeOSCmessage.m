function buff = writeOSCmessage(msg, dat)
%WRITEOSCMESSAGE Encode an OSC message.
%   B = WRITEOSCMESSAGE(MSG, DAT) encodes the OSC message (MSG) and data
%   (DAT) into a character buffer (B) for sending across UDP.
%
%   See also SENDOSCMESSAGE, READOSCMESSAGE.

%   ==============================================================================
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
%   ==============================================================================

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