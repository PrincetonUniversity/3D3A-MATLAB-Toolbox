function erb = fc2erb(fc,n)
%fc2erb Equivalent rectangular bandwidth (ERB) at a center frequency.
%   ERB = fc2erb(Fc) computes the ERB at a center frequency Fc, given in
%   Hz.
%
%   ERB = fc2erb(Fc,n) uses the nth order polynomial approximation given by
%   Moore and Glasberg. Accepts n = 1 or n = 2 only.

if nargin < 2
    n = 1;
end

fc = fc/1000; % convert Hz to kHz

switch n
    case 1
        % The approximation is applicable at moderate sound levels and for
        % values of fc between 0.1 and 10 kHz.
        erb = 24.7*(4.37*fc + 1);
    case 2
        % The approximation is based on the results of a number of
        % published simultaneous masking experiments and is valid from 0.1
        % to 6.5 kHz.
        erb = 6.23*(fc.^2) + 93.39*fc + 28.52;
    otherwise
        error('No polynomial approximation is known for n = %g',n)
end