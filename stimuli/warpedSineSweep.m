function x = warpedSineSweep(Xmag, Fs, T, varargin)
%warpedSineSweep Generate warped sine sweep.
%   x = warpedSineSweep(Xmag, Fs, T) generates a warped SS at sampling rate
%   Fs that has a magnitude spectrum Xmag and duration T seconds.
%
%   Beyond the 3 required inputs, warpedSineSweep accepts the following
%   optional argument(s):
%       {'pad', Tpad} - pads the sweep with Tpad seconds of zeros before
%           and after the sweep.
%
%       {'pad', [Tpre, Tpost]} - pads the sweep with Tpre seconds of zeros
%           before and Tpost seconds of zeros after the sweep.
%   
%   Equations from Ochiai and Kaneda (2013) "A Recursive Adaptive Method of
%   Impulse Response Measurement with Constant SNR over Target Frequency
%   Band," JAES.

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