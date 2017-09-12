function f = k2f(k)
%k2f Convert angular wavenumber to temporal frequency in Hz.
%   f = k2f(k) returns temporal frequency f given by k*C/(2*PI), where C is
%   the speed of sound.
%
%   See also f2k.

f = k*getSoundSpeed()/(2*pi);

end