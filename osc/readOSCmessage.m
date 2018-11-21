function [addr, args] = readOSCmessage(buff)
%READOSCMESSAGE Parse an OSC message.
%   [MSG, DAT] = READOSCMESSAGE(B) parses the OSC character buffer (B) into
%   the OSC message (MSG) and the corresponding data (DAT).
%
%   See also WRITEOSCMESSAGE.

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