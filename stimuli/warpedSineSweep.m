function x = warpedSineSweep(Xmag, Fs, T, varargin)
%WARPEDSINESWEEP Generate a warped sine sweep.
%   X = WARPEDSINESWEEP(XMAG, FS, T) generates a warped sine sweep at
%   sampling rate FS that has a magnitude spectrum XMAG and duration T
%   seconds.
%
%   WARPEDSINESWEEP(XMAG,FS,T,Name1,Value1,Name2,Value2, ...) 
%   Specifies optional comma-separated pairs of Name,Value arguments, 
%   where Name is the argument name and Value is the corresponding value. 
%   Name must appear inside single quotes (' '). You can specify several 
%   name and value pair arguments in any order as Name1,Value1,...,NameN,
%   ValueN.  Valid Name,Value arguments are as follows:
%
%   'pad'   Duration (in seconds) of zeros before and after the sweep. If
%           given as a scalar, TPAD, this option pads the sweep with TPAD
%           seconds of zeros before and after the sweep. If given as a
%           pair, [TPRE, TPOST], this option pads the sweep with TPRE
%           seconds of zeros before and TPOST seconds of zeros after the
%           sweep.

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

%   References:
%     [1] Ochiai and Kaneda (2013) A Recursive Adaptive Method of Impulse
%         Response Measurement with Constant SNR over Target Frequency Band.

narginchk(3,inf);

% Check for padding option
indx = find(strcmpi(varargin,'pad'),1,'first');
if isempty(indx)
    Tpre = 0;
    Tpost = 0;
else
    Tpad = varargin{indx+1};
    switch numel(Tpad)
        case 1
            Tpre = Tpad;
            Tpost = Tpad;
        case 2
            Tpre = Tpad(1);
            Tpost = Tpad(2);
    end
end

% Compute sweep
L = Fs*T;
Xmag = abs(fft(ifft(Xmag,'symmetric'),L));

a1 = L/(sum(Xmag(1:(1+L/2)).^2) - Xmag(1).^2); % Eq. (12)
D = a1*(cumsum(Xmag(1:(1+L/2)).^2) - Xmag(1).^2); % Eq. (10)
Xphase = (2*pi/L)*cumsum(D); % Eq. (13)

X = zeros(L,1);
X(1:(1+L/2)) = Xmag(1:(1+L/2)).*exp(-1i*Xphase); % Eq. (14)
sweep = ifft(X,L,'symmetric');

x = [zeros(Fs*Tpre,1); sweep/max(abs(sweep)); zeros(Fs*Tpost,1)];

end