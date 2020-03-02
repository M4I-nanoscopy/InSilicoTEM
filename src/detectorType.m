function [mtf, dqe, cf, readn, darkn] = detectorType(params,convfact)
%detectorType: Reads in the paratemres of already characterized detectors
% SYNOPSIS:
% [mtf, dqe, cf, readn, darkn] = detectorType(params)
%
% PARAMETERS:
%   params: structure containing various input parameters (e.g. camera type, voltage, directory)
%
% OUTPUT:
%    mtf: modulation transfer function
%    dqe: detective quantum efficieny
%  readn: readout noise
%  darkn: dark current noise

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic
cam          = params.cam.type;
Voltage      = params.acquis.Voltage;
camdir       = params.camdir;
     
s0=  'MTF is not available for this voltage: %03d keV. We will inter(extra)polate from the MTFs measured at %03d and %03d keV\n';
s01= 'MTF is not available for this voltage: %03d keV. We will extrapolate from the MTF measured at %03d keV\n';
s02= 'MTF is not available for this voltage: %03d keV. We will inter(extra)polate from the MTFs measured at %03d, %03d and %03d keV\n';

if (strcmp(cam,'Eagle4k') && Voltage/1000~=120) 
    s1=sprintf(s01, Voltage/1000, 120);
    warning(s1);
elseif (strcmp(cam,'US4000') && (Voltage/1000~=120 && Voltage/1000~=200))
    s1=sprintf(s0, Voltage/1000, 120, 200);
    warning(s1);  
elseif strcmp(cam,'FalconI') && (Voltage/1000~=200 && Voltage/1000~=300)
    s1=sprintf(s0, Voltage/1000, 200, 300);
    warning(s1);
elseif strcmp(cam,'FalconIII_Linear') && (Voltage/1000~=200 && Voltage/1000~=300)
    s1=sprintf(s0, Voltage/1000, 200, 300);
    warning(s1);
elseif strcmp(cam,'FalconIII_EC') && (Voltage/1000~=200 && Voltage/1000~=300)
    s1=sprintf(s0, Voltage/1000, 200, 300);
    warning(s1);
elseif strcmp(cam,'US1000GIF') && (Voltage/1000~=80 && Voltage/1000~=200 && Voltage/1000~=300)
    s1=sprintf(s02, Voltage/1000, 80, 200, 300);
    warning(s1);
elseif strcmp(cam,'custom')
    warning('This camera has not been characterized. Please provide a MTF and DQE. Or switch camera type to ''custom''');     
end
      
if  strcmp(cam,'ideal')
        cf    = 1000;    
        mtf   = dip_image(ones(256,256));
        dqe  = mtf^2;
        readn = 0;
        darkn = 0;
            
elseif strcmp(cam,'perfect') 
        cf    = 1000;
        mtf   = sinc(xx/size(xx,1))*sinc(yy/size(yy,1));
        dqe  = mtf^2;
        readn = 0;
        darkn = 0;   
            
elseif strcmp(cam,'US1000GIF') % Gatan US1000. Data available for 80, 200 and 300 kV. For the other voltages interpolation is used                 
       if   Voltage == 80e3 
            cf     = 28;
            bz     = load([camdir filesep 'MTF_US1000GIF_080']);
            mtf    = bz.mtf;
            if exist([camdir filesep 'DQE_US1000GIF_080.mat'])
                npsz   = load([camdir filesep 'DQE_US1000GIF_080']);
                dqe   = npsz.dqe;
            else
               warning('DQE not available. Use only MTF.')
               dqe = mtf^2;
            end
            readn  = 3;
            darkn  = 0.11;            
       elseif Voltage == 200e3 
            cf     = 15; %check
            bz     = load([camdir filesep 'MTF_US1000GIF_200']);
            mtf    = bz.mtf;
            if exist([camdir filesep 'DQE_US1000GIF_200.mat'])
                npsz   = load([camdir filesep 'DQE_US1000GIF_200']);
                dqe   = npsz.dqe;
            else
                warning('DQE not available. Use only MTF.')
                dqe = mtf^2;
            end
            readn  = 3;
            darkn  = 0.11; 
       elseif Voltage == 300e3 
            cf     = 6.6;
            bz     = load([camdir filesep 'MTF_US1000GIF_300']);
            mtf    = bz.mtf;
            if exist([camdir filesep 'DQE_US1000GIF_300.mat'])
                npsz   = load([camdir filesep 'DQE_US1000GIF_300']);
                dqe   = npsz.dqe;
            else
                warning('DQE not available. Use only MTF.')
                dqe = mtf^2;
            end
            readn  = 3;
            darkn  = 0.11; 
       else19
            cf     = (sqrt(80)*28*sqrt(200)*15*sqrt(300)*6.6)^(1/3)/sqrt(Voltage/1000); 
            bz     = load([camdir filesep 'MTF_US1000GIF_080']);
            mtf80  = bz.mtf;
            bz     = load([camdir filesep 'MTF_US1000GIF_200']);
            mtf200 = bz.mtf;
            bz     = load([camdir filesep 'MTF_US1000GIF_300']);
            mtf300 = bz.mtf;
            mtf    = (sqrt(80)*mtf80*sqrt(200)*mtf200*sqrt(300)*mtf300)^(1/3)/sqrt(Voltage/1000);
            mtf    = mtf/max(mtf);
            if exist([camdir filesep 'DQE_US1000GIF_080.mat'])&& exist([camdir filesep 'DQE_US1000GIF_200.mat']) && exist([camdir filesep 'DQE_US1000GIF_300.mat'])
                npsz   = load([camdir filesep 'DQE_US1000GIF_080']);
                dqe80 = npsz.dqe;
                npsz   = load([camdir filesep 'DQE_US1000GIF_200']);
                dqe200= npsz.dqe;
                npsz   = load([camdir filesep 'DQE_US1000GIF_300']);
                dqe300= npsz.dqe; 
                dqe   = (sqrt(80)*dqe80*sqrt(200)*dqe200*sqrt(300)*dqe300)^(1/3)/sqrt(Voltage/1000);
            else
                warning('Noise power spectra not available. Use only MTF.')
                dqe = mtf^2;
            end
            readn  = 3;
            darkn  = 0.11; 
       end  
       
elseif strcmp(cam,'Eagle4k')    
       if Voltage == 120e3 
            cf     = 100; %check 
            bz     = load([camdir filesep 'MTF_Eagle_120']);
            mtf    = bz.mtf;
            if exist([camdir filesep 'DQE_Eagle_120.mat'])
                npsz   = load([camdir filesep 'DQE_Eagle_120']);
                dqe   = npsz.dqe;
            else
                warning('DQE not available. Use only MTF.')
                dqe = mtf^2;
            end
            readn  = 7;
            darkn  = 1.14; 
       else
            cf     = sqrt(120)*100/sqrt(Voltage/1000);
            bz     = load([camdir filesep 'MTF_Eagle_120']);
            mtf120 = bz.mtf;
            mtf    = sqrt(120)*mtf120/sqrt(Voltage/1000);
            mtf    = mtf/max(mtf);
            if exist([camdir filesep 'DQE_Eagle_120.mat'])
                npsz   = load([camdir filesep 'DQE_Eagle_120']);
                dqe120= npsz.dqe;
                dqe    = sqrt(120)*dqe120/sqrt(Voltage/1000);
            else
                warning('DQE not available. Use only MTF.')
                dqe = mtf^2;
            end           
            readn  = 7;
            darkn  = 1.14; 
       end  
       
elseif strcmp(cam,'US4000')    
       if  Voltage == 120e3 
            cf      = 34;
            bz      = load([camdir filesep 'MTF_US4000_120']);
            mtf     = bz.mtf;
            if exist([camdir filesep 'DQE_US4000_120.mat'])
                npsz    = load([camdir filesep 'DQE_US4000_120']);
                dqe    = npsz.dqe;
            else
                warning('DQE not available. Use only MTF.')
                dqe = mtf^2;
            end
            readn   = 3;
            darkn   = 0.11; 
      elseif Voltage == 200e3 
            cf       = 23;
            bz       = load([camdir filesep 'MTF_US4000_200']);
            mtf      = bz.mtf;
            if exist([camdir filesep 'DQE_US4000_200.mat'])
                npsz   = load([camdir filesep 'DQE_US4000_200']);
                dqe    = npsz.dqe;
            else
                warning('DQE not available. Use only MTF.')
                dqe = mtf^2;
            end
            readn    = 3;
            darkn    = 0.11;             
      else
            cf       = sqrt(sqrt(120)*34*sqrt(200)*23)/sqrt(Voltage/1000);
            bz       = load([camdir filesep 'MTF_US4000_120']);
            mtf120   = bz.mtf;
            bz       = load([camdir filesep 'MTF_US4000_200']);
            mtf200   = bz.mtf;
            mtf      = sqrt(sqrt(120)*mtf120*sqrt(200)*mtf200)/sqrt(Voltage/1000);
            mtf      = mtf/max(mtf);
            if exist([camdir filesep 'DQE_US4000_120.mat'])&& exist([camdir filesep 'DQE_US4000_200.mat'])
                npsz     = load([camdir filesep 'DQE_US4000_120']);
                dqe120  = npsz.dqe;
                npsz     = load([camdir filesep 'DQE_US4000_200']);
                dqe200  = npsz.dqe;
                dqe     = sqrt(sqrt(120)*dqe120*sqrt(200)*dqe200)/sqrt(Voltage/1000);
            else
                warning('Noise power spectra not available. Use only MTF.')
                dqe = mtf^2;
            end            
            readn    = 3;
            darkn    = 0.11; 
       end  
elseif strcmp(cam,'FalconIII_Linear')
       if Voltage  == 200e3
            cf      = 19;
            bz      = load([camdir filesep 'MTF_FalconIII_Linear_200']);
            mtf     = bz.mt19f;
            npsz    = load([camdir filesep 'DQE_FalconIII_Linear_200']);
            dqe     = npsz.dqe;
            readn   = 3;
            darkn   = 1.6;
       elseif Voltage == 300e3
            cf      = 19;
            bz      = load([camdir filesep 'MTF_FalconIII_Linear_300']);
            mtf     = bz.mtf;
            npsz    = load([camdir filesep 'DQE_FalconIII_Linear_300']);
            dqe     = npsz.dqe;
            readn   = 3;
            darkn   = 1.6;
       end
       
elseif strcmp(cam,'FalconIII_EC')
       if Voltage  == 200e3
            cf      = convfact; % Conversion factor ?
            bz      = load([camdir filesep 'MTF_FalconIII_EC_200']);
            mtf     = bz.mtf;
            if  exist([camdir filesep 'DQE_FalconIII_EC_200.mat'])
                npsz   = load([camdir filesep 'DQE_FalconIII_EC_200']);
                dqe   = npsz.dqe;
            else
                warning('DQE not available. Use only MTF.')
                dqe = mtf^2;
            end
            readn   = 3;
            darkn   = 1.6; 
       elseif Voltage == 300e3
            cf      = convfact; % Conversion factor
            bz      = load([camdir filesep 'MTF_FalconIII_EC_300']);
            mtf     = bz.mtf;
            npsz    = load([camdir filesep 'DQE_FalconIII_EC_300']);
            dqe     = npsz.dqe;
            readn   = 3;
            darkn   = 1.6;
       end
       
elseif strcmp(cam,'FalconI')    
       if  Voltage == 200e3 
            cf      = 120;
            bz      = load([camdir filesep 'MTF_Falcon_200']);
            mtf     = bz.mtf;
            npsz    = load([camdir filesep 'DQE_Falcon_200']);
            dqe     = npsz.dqe;
            readn   = 3;
            darkn   = 0.11; 
      elseif Voltage == 300e3 
            cf       = 100;
            bz       = load([camdir filesep 'MTF_Falcon_300']);
            mtf      = bz.mtf;
            npsz     = load([camdir filesep 'DQE_Falcon_300']);
            dqe      = npsz.dqe;
            readn    = 2;
            darkn    = 0.11;             
      else
            cf       = 100;
            bz       = load([camdir filesep 'MTF_Falcon_200']);
            mtf200   = bz.mtf;
            bz       = load([camdir filesep 'MTF_Falcon_300']);
            mtf300   = bz.mtf;
            mtf      = sqrt(sqrt(200)*mtf200*sqrt(300)*mtf300)/sqrt(Voltage/1000);
            mtf      = mtf/max(mtf);
            if exist([camdir filesep 'DQE_Falcon_200.mat'])&& exist([camdir filesep 'DQE_Falcon_300.mat'])
                npsz     = load([camdir filesep 'DQE_Falcon_200']);
                dqe200  = npsz.dqe;
                npsz     = load([camdir filesep 'DQE_Falcon_300']);
                dqe300  = npsz.dqe;
                dqe     = sqrt(sqrt(200)*dqe200*sqrt(300)*dqe300)/sqrt(Voltage/1000);
            else
                warning('Noise power spectra not available. Use only MTF.')
                dqe = mtf^2;
            end 
            readn    = 3;
            darkn    = 0.11; 
       end   
elseif strcmp(cam,'custom') 
    warning('Parameters of this camera are not known. Please characterize it with the supplementary software or simulate it as EMG (exponentially modified Gaussian) by setting the flag "params.cam.GenMTFasEMG" ');
    if params.cam.GenMTFasEMG
        mtf = SimMTFasEMG(params);
        if params.cam.DQEflag           
            dqe = params.cam.dqe0*(mtf/mtf^params.cam.mtf_exp4nnps+eps)^2;
            cf = 1000*sqrt(params.cam.dqe0./(params.acquis.dose_on_sample*(params.acquis.pixsize*1e10)^2));
        else
            dqe = 1+0*mtf;
            cf = 1000*sqrt(1./(params.acquis.dose_on_sample*(params.acquis.pixsize*1e10)^2));
        end        
    end
    readn = params.cam.readn;
    darkn = params.cam.darkn; 
else
    error('This camera type is not known\n')
end
    % Avoid possible NaN
            mtf(mtf==0) = 2;
%             dipshow(mtf)
            minmtf=min(mtf);
%             dipshow(minmtf)
            mtf= mtf*(mtf~=2)+dip_image(minmtf*double((mtf==2)));
            dqe(dqe==0) = 2;
            mindqe=min(dqe);
            dqe= dqe*(dqe~=2)+dip_image(mindqe*double((dqe==2)));
            
if isfield(params.cam, 'nnps_empirical') && ~strcmp(cam,'ideal') && ~strcmp(cam,'perfect')
    if params.cam.nnps_empirical
        dqe  = params.cam.dqe0*(mtf/mtf^params.cam.mtf_exp4nnps+eps)^2;
    else
        dqe  = params.cam.dqe0*(1+0*mtf);
    end
    if strcmp(cam,'custom')
       cf = 1000*sqrt(params.cam.dqe0./(params.acquis.dose_on_sample*(params.acquis.pixsize*1e10)^2));
    end
end

% MTF
mtf    = resample(mtf,params.proc.N/256);
mtf1   = cut(mtf,params.proc.N/params.cam.bin); 
mtf    = resample(mtf1, params.cam.bin); 

% DQE
dqe    = resample(dqe,params.proc.N/256); 
dqe1   = cut(dqe,params.proc.N/params.cam.bin);
dqe    = resample(dqe1, params.cam.bin); 
