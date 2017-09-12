function lambda = f2lambda(f)
%f2lambda Convert temporal frequency in Hz to wavelength.
%   lambda = f2lambda(f) returns wavelength lambda given by C/f, where C is
%   the speed of sound.

lambda = getSoundSpeed()./f;

end