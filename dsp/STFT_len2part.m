function nparts = STFT_len2part(len, partLen, novlp, padFlag)
%STFT_LEN2PART Number of partitions needed for a signal.
%   NPARTS = STFT_LEN2PART(LEN,PARTLEN,NOVERLAP) returns the number of
%   partitions NPARTS, each of length PARTLEN and overlapping by NOVERLAP
%   samples, needed to partition a signal of length LEN.
%
%   NPARTS = STFT_LEN2PART(LEN,PARTLEN,NOVERLAP,PAD) optionally assumes the
%   signal will be padded to at least LEN+PARTLEN and rounds NPARTS up if
%   PAD evaluates to true. By default, the signal is not assumed to be
%   padded and NPARTS is rounded down.
%
%   See also STFT_PART2LEN, GETFORWARDSTFT, GETINVERSESTFT.

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

if nargin < 4 || isempty(padFlag)
    padFlag = false;
end

hop = partLen - novlp;
if padFlag
    nparts = ceil((len + hop) / hop);
else
    nparts = fix((len - novlp) / hop);
end

end