function abl = motionBlur(a,params2) 
%motionBlur introduces extra blurring of the potential which can be ralted
% to the beam-indeuced motion
% SYNOPSIS:
% abl = motionBlur(a,params2) 
%
% PARAMETERS:
%        a: input potental of the particle (without motion assumed) 
%  params2: structure containing various input paramters (e.g. motion, voxel size)
%
% OUTPUT:
%  PartPot: Interaction potential of the particle 
%  pout   : possibly changed parameters (number of particles)

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

    a = dip_image(a);
  aft = ft(a);
 Btot = (8*pi^2*(params2.spec.motblur)^2)/4; %devided by 4 in order to bring ksi (frequency) from crystalography to q from EM. 
  abl = real(ift(aft*exp(-Btot*(rr(aft,'freq')/params2.spec.voxsize).^2)));
  abl = double(abl);