
function lambda = wavlen(V)
% relativistic wavelength gievn the accelration voltage V
h = 6.626068e-34;
e = 1.60217646e-19;
m = 9.10938188e-31;
c = 2.99792458e8;

lambda = h/sqrt(e*V*m*(e/m*V/c^2 + 2 ));

end