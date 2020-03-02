function psi_multislice  = multislice(psi_t,Nm,n,lambda,q_m,dzprop)
%multislice Performs multislice calculations 
% SYNOPSIS:
% psi_multislice  = multislice(psi_t,Nm,n,lambda,q_m,dzprop)
%
% PARAMETERS:
%   psi_t: the phase grating 
%      Nm: the sieze of the image
%       n: number of slices
%  lambda: wavelength of the incoming electrons
%      qm: frequencies in the Fourier domain
%  dzprop: propagation distance of the Fresnel propagator
%
% OUTPUT:
%  psi_multislice: electron wave exiting the specimen

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

% Fresnel propagator
P = exp(-1i*pi*lambda*(q_m.^2)*dzprop); % Fresnel propagator

% MULTISLICE
psi_multislice = newim(Nm,Nm,'dcomplex')+1;
for ii = 1:n
    psi_multislice = ift(ft(psi_multislice*squeeze(psi_t(:,:,ii-1)))*P);
    %psi_multislice = ift(ft(psi_multislice*squeeze(psi_t(:,:,ii-1)))*P(dzprop));
end