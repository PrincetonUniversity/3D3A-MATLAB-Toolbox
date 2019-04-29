function x = exponentialSineSweep(f1, f2, Fs, T, varargin)
%EXPONENTIALSINESWEEP Generate an exponential sine sweep (ESS).
%   X = EXPONENTIALSINESWEEP(F1, F2, FS, T) generates an ESS at sampling
%   rate FS that sweeps from F1 to F2 in T seconds.
%
%   EXPONENTIALSINESWEEP(F1,F2,FS,T,Name1,Value1,Name2,Value2, ...) 
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
%
%   'type'  Type of ESS to generate. Supported types are 'standard' (Farina
%           ESS) and 'phase', which generates a phase-controlled ESS. Note
%           that the specified start frequency F1 and the duration T are
%           only satisfied approximately for 'phase' type sweeps.

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
%     [1] Farina (2000) Simultaneous Measurement of Impulse Response and
%         Distortion with a Swept-Sine Technique.
%     [2] Vetter and di Rosario (2011) ExpoChirpToolbox: a Pure Data
%         implementation of ESS impulse response measurement.

narginchk(4,inf);

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

% Check for sweep method
indx = find(strcmpi(varargin,'type'),1,'first');
if isempty(indx)
    TYPE = 'standard';
else
    TYPE = varargin{indx+1};
end

% Compute sweep
switch lower(TYPE)
    case {'standard','farina'}
        w1 = 2*pi*f1/Fs;
        w2 = 2*pi*f2/Fs;
        N = Fs*T;
        
        k = (0:(N-1)).';
        sweep = sin((w1*N/log(w2/w1))*((w2/w1).^(k/N)-1));
    case {'phase','pc','phase controlled','phase-controlled'}
        P = round(log2(f2/f1));
        w2 = 2*pi*f2/Fs;
        w1 = w2/(2^P);
        
        L = 2*pi*log(2^P)*round(w1*T*Fs/(2*pi*log(2^P)))/w1;
        N = round(L);
        
        k = (0:(N-1)).';
        sweep = sin((w1*L/log(2^P))*((2^P).^(k/N)));
    otherwise
        error('No such sweep method is known.');
end

x = [zeros(Fs*Tpre,1); sweep; zeros(Fs*Tpost,1)];

end
        