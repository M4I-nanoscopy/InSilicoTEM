function Lambda_in = in_MFP(params, rho, Mw)
%in_MFP Calculates the inelastic mean free path given the composition of
%the specimen and the accelarating voltage according to the formula 8 in
%the article
% SYNOPSIS:
% Lambda_in = in_MFP(params, rho, Mw)
%
% PARAMETERS:
%   params: structure containing various input paramters (e.g. Voltage)
%
% OUTPUT:
%  Lambda_in: inelastic mean free path

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic
nc = phys_const();
U0 = params.acquis.Voltage; % incident energy
E0 = nc.me*nc.c^2/nc.el;    %rest mass energy
E_loss = 20;                % eV mean plasmon losses
beta2 = 1-(E0/(U0+E0))^2;
beta2_100 = 1-(E0/(100e3+E0))^2;

if     Mw == 18 % water
    ZO=8;
    sigma_inH=8.8e-6*beta2_100.*log(beta2.*(U0+E0)/(E_loss/2))/(beta2.*log(beta2_100.*(100e3+E0)/(E_loss/2)));
    sigma_inO=1.5*1e-6*ZO^0.5./beta2.*log(beta2.*(U0+E0)/(E_loss/2));
    sigma_in=2*sigma_inH+sigma_inO;
    
elseif Mw == 12 %carbon
    ZC=6;
    sigma_in=1.5*1e-6*ZC^0.5./beta2.*log(beta2.*(U0+E0)/(E_loss/2));
    
elseif Mw == 7.2 % protein   
    sigma_in=0.82*1e-4*beta2_100.*log(beta2.*(U0+E0)/(E_loss/2))/(beta2.*log(beta2_100.*(100e3+E0)/(E_loss/2)));
end

Nconc   = rho*1000*nc.Na/(Mw/1000);
Lambda_in=1./(Nconc*sigma_in*1e-18);


