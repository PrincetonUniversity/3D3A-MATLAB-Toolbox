function x = exponentialSineSweep(f1, f2, Fs, T, varargin)
%exponentialSineSweep Generate exponential sine sweep (ESS).
%   x = exponentialSineSweep(f1, f2, Fs, T) generates an ESS at sampling
%   rate Fs that sweeps from f1 to f2 in T seconds.
%
%   Beyond the 4 required inputs, exponentialSineSweep accepts the
%   following optional argument(s):
%       {'pad', Tpad} - pads the sweep with Tpad seconds of zeros before
%           and after the sweep.
%
%       {'pad', [Tpre, Tpost]} - pads the sweep with Tpre seconds of zeros
%           before and Tpost seconds of zeros after the sweep.
%
%       {'method', 'phase'} - generates a phase-controlled ESS. Note that
%           the specified start frequency f1 and the duration T are only
%           satisfied approximately.

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
indx = find(strcmpi(varargin,'method'),1,'first');
if isempty(indx)
    METHOD = 'standard';
else
    METHOD = varargin{indx+1};
end

% Compute sweep
switch lower(METHOD)
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
        