function params2 = parsePar(params)
%parsePar Parses input parameters.
% SYNOPSIS:
% [params2] = parsePar(params)
%
% PARAMETERS:
%   params: structure containing various input physical and processing
%   parameters (fields: proc, spec, inter, mic, acuis, cam, disp)
%
% OUTPUT:
%  params2: structure containing various output physical and processing parameters

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

params2 = params;
  
    
% ------------------------- Detector parameters-------
if ispc; sep='\'; else sep='/';end
params.camdir = [pwd sep 'MTFs']; % place the estiamted MTFs of the detector into MTFs subfolder

 % In case params.cam.type = 'custom' and params.cam.DQEflag = 1; specify the parameters of the dqe
if  strcmp(params.cam.type,'custom')
    params.cam.readn          =   3; % Readout noise in [ADUs]
    params.cam.darkn          = 0.7; % Dark current noise in [ADUs]
    if params.cam.DQEflag ~= 0
         params.cam.nnps_empirical = 1;   % flag: to simulate nnps (normalized noise power spectrum) as mtf^exponent in case nnps is not available
         params.cam.mtf_exp4nnps   = 3/8; % in case dqe is not available define empiricly how much better nnps should be compared to the mtf 
         params.cam.dqe0           = 0.4; % DQE at zero frequency
    end
end
% load or generate mtf and dqe
[mtf, dqe, cf, readn, darkn] = detectorType(params);
    params2.cam.mtf   = mtf;
    params2.cam.dqe   = dqe;
    params2.cam.cf    = cf; 
    params2.cam.readn = readn;
    params2.cam.darkn = darkn; 


%------------------------- Series (possible only one at the time)
params2.seriesout='none';
% Note: At this moment it is possible to generate defocus, flux or tilt series (but not simultaneously) 
if length(params.acquis.tilt)>1
    params2.seriesout = 'tilt';
    if length(params.acquis.dose_on_sample)>1 || length(params.acquis.df)>1
    warning('Ignoring other entries than the first for defocus and/or flux. It is possible to generate only one series at the time. Generating a tilt series...')  
    params2.acquis.df = params.acquis.df(1);
    params2.acquis.dose_on_sample = params.acquis.dose_on_sample(1);
    end
elseif length(params.acquis.df)>1
    params2.seriesout ='defocus';
    if length(params.acquis.dose_on_sample)>1 
    warning('Ignoring other entries than the first for flux. It is possible to generate only one series at the time. Generating a defocus series...')  
    params.acquis.dose_on_sample = params.acquis.dose_on_sample(1);
    end
elseif length(params.acquis.dose_on_sample)>1
    params2.seriesout = 'dose'; 
end

% --------------------- Other paramters    
nc = phys_const();
params2.inter.lambda       = wavlen(params.acquis.Voltage); % wavelenght
           relmass         = nc.me + nc.el*params.acquis.Voltage/(nc.c^2); % relativistic mass
params2.inter.sig_transfer = 2*pi*relmass*nc.el*params2.inter.lambda/(nc.h^2); % interaction constant
params2.acquis.dose        = params.acquis.dose_on_sample*(params.acquis.pixsize*1e10)^2; %el/pixel
params2.spec.Vwat          = 4.5301; % The potential value (constant) of the amourphous ice
params2.proc.N             = params.proc.N/params.cam.bin; 
params2.spec.voxsize       = 1; %% Voxel size of the potential map in [A](Default: 1A, the accuracy decreases if voxsize too large)

% Inelastics
params2.spec.imagpot2specm   = 'amorIce'; % Options: 'amorC' and 'amorIce', 'graphite' %used when params.amplconsem =2
params2.spec.imagpot1Q       = 0.075; % Required if params.amplconsem = 1
if     params.spec.imagpot == 1
          	params2.spec.potenampl = params2.spec.imagpot1Q;
elseif params.spec.imagpot == 2 || params.spec.imagpot == 3
    switch params2.spec.imagpot2specm
        case 'amorC'
           Lambda_in = in_MFP(params, 1.8, 12); %inelastic mean free path
        case 'graphite'
           Lambda_in = in_MFP(params, 2.3, 12);
        case 'amorIce'
           Lambda_in = in_MFP(params, 0.93, 18);
    end 
  params2.spec.potenampl = 1/(2*params2.inter.sig_transfer*Lambda_in); 
    if params.spec.imagpot == 3
       Lambda_in_prot = in_MFP(params, 1.35, 7.2);
       params2.spec.proteinampl=1/(2*params2.inter.sig_transfer*Lambda_in_prot);
    end
end
    
    
 
    
    

   