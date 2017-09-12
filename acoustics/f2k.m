function k = f2k(f)
%f2k Convert temporal frequency in Hz to angular wavenumber.
%   k = f2k(f) returns angular wavenumber k given by 2*PI*f/C, where C is
%   the speed of sound.
%
%   See also f2k.

k = 2*pi*f/getSoundSpeed();

end