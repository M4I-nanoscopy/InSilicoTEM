function imfinal = DetectSim(btot_i, params2)
%DetectSim simulates the effects of the detector (mtf, dqe, noise)
%  
% SYNOPSIS:
% imfinal = DetectSim(btot_i, params2)
%
% PARAMETERS:
%  btot_i: intensity in the image just before hitting the detector (noiseless) 
% params2: structure containing various input physical and processing parameters
%
% OUTPUT:
% imfinal: final noisy image after the detector

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

b = ft(btot_i);

% MTF 
mtffin = params2.cam.mtf;

% DQE
dqefin = params2.cam.dqe;
    
if params2.cam.DQEflag == 1 
    nnps = mtffin^2/dqefin;
    nnps(nnps==0) = max(min(nnps),1e-7);
    Sn = b*mtffin./sqrt(nnps);
    btot_i2 = real(ift(Sn)); 
    ImPoission = dip_image(noise(double(btot_i2*params2.influx),'poisson',1)); %poisson noise
    d2f = params2.cam.cf*ft(ImPoission)*sqrt(nnps); 
else     % only MTF
    d1  = dip_image(noise(double(btot_i*params2.influx),'poisson',1)); %poisson noise
    d2f = params2.cam.cf*squeeze(ft(d1)).*mtffin;
end

%readout and dark noise
 readsim = noise((params2.cam.readn),'gaussian'); %readout noise has a gaussian distribution 
 darksim = 0; %dark current noise has a poisson distribution, usually an order of magnitude smaller than readout noise
 imfinal = real(ift(squeeze(d2f))) + readsim + darksim; % the final single tilt image

% display mtf or dqe
if params2.disp.mtfdqe && params2.cam.DQEflag == 1 
     dqe = 1*mtffin.^2/nnps;
     sz2 = size(dqe,1)/2;
     x   = [1:sz2]./sz2;
     figure(40)
     plot(x, mtffin(sz2:sz2*2-1, sz2), 'b-', x, dqe(sz2:sz2*2-1, sz2), 'r-.')
     legend('MTF', 'DQE')
     xlabel('Nyquist')
     xlim([ 0 1.1 ])
     ylim([ 0 1.1 ])
elseif params2.disp.mtfdqe && params2.cam.DQEflag == 0 
     nnps = mtffin^2;
     dqe  = 1*mtffin.^2/nnps;
     sz2  = size(mtffin,1)/2;
     x    = [1:sz2]./sz2;
     figure(40)
     plot(x, mtffin(sz2:sz2*2-1, sz2), 'b-', x, dqe(sz2:sz2*2-1, sz2), 'r-.')
     legend('MTF', 'DQE')
     xlabel('Nyquist')
     xlim([ 0 1.1 ])
     ylim([ 0 1.1 ])
end
        
        