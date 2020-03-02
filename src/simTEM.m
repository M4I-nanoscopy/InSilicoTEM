function imStructOut= simTEM(InputVol, params2)
%simTEM simulates the full image formation given the input interaction
% potential and acquisition parameters (dose, optics and detector) 
%
% SYNOPSIS:
% [imStructOut,paramsout]= simTEM(InputVol, params2)
%
% PARAMETERS:
%  InputVol: Interaction potenital specimen volume 
%   params2: Structure containing various input simulation paramters
%
% OUTPUT:
% imStructOut: Structure containg output images (or stacks): noisy, noiseless and exit wave image 
%              Optional: ctf and mtf (if params.disp.ctf or params.disp.mtfdqe are set)

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic
%
%
voxSz = params2.acquis.pixsize;% the voxel size
% prealocate memory for the stack (series)
nTiltAngles = length(params2.acquis.tilt);

if strcmp(params2.seriesout,'tilt')
    Nseries= nTiltAngles;
elseif strcmp(params2.seriesout,'defocus')
    Nseries= length(params2.acquis.df);
elseif strcmp(params2.seriesout,'dose')
    Nseries= length(params2.acquis.dose_on_sample);
else 
    Nseries=1;
end
series        = newim(params2.proc.N,params2.proc.N,Nseries); 
btot_i        = newim(params2.proc.N,params2.proc.N,Nseries);
extprojstack  = newim(params2.proc.N,params2.proc.N,Nseries, 'complex');

poten0 = InputVol; clear InputVol;
thickness = voxSz*size(poten0,3);

switch params2.inter.type
    case{'pa+wpoa', 'pa', 'wpoa', 'tpga'}
        if params2.proc.cores>1
           parpool(params2.proc.cores);
        end
        poten0 = permute(poten0,[1 3 2]);
    % calculate the projected potential (neglects the thickness of the specimen)
    switch params2.inter.type
        case {'pa+wpoa','pa'}
            poten0_stack = FreqDomRadon_2d(poten0,params2.acquis.tilt);
    % in case of wpoa and tpga the Ewald sphere (parabola) is sampled and propagation through the specimen thickness is taken into account 
        case {'wpoa','tpga'}
            poten0_stack = FreqDomRadon_Ewald_2d(poten0,params2.acquis.tilt, params2);
    end
    if params2.proc.cores>1
        ppc = gcp;
        delete(ppc);
    end
    poten0_stack = permute(poten0_stack,[1 3 2]);
    poten0_stack = cut(poten0_stack,[size(InputVol,1) size(InputVol,2) size(poten0_stack,3)]);
    
    % imaginary part of the potential (amplitude contrast)
    if params2.spec.imagpot == 1
       poten_stack = poten0_stack + 1i*params2.spec.potenampl*poten0_stack;
    elseif params2.spec.imagpot == 2 || params2.spec.imagpot == 3
       imagPotProj = double((newim(poten0_stack)+1)*dip_image(reshape(params2.spec.potenampl*round(thickness/(voxSz))./cos(params2.acquis.tilt), [1 1 nTiltAngles])));
       poten_stack = poten0_stack + 1i*imagPotProj;
    else
       poten_stack = poten0_stack;
    end
    
    % transmission function
    if strcmp(params2.inter.type, 'pa')|| strcmp(params2.inter.type, 'tpga') 
        psi_exit   = exp(1i*params2.inter.sig_transfer*poten_stack*voxSz);
    % the first term (corresponds to the single scattering)    
    elseif strcmp(params2.inter.type, 'pa+wpoa')|| strcmp(params2.inter.type, 'wpoa')
        psi_exit   = 1+1i*params2.inter.sig_transfer*poten_stack*voxSz;
    end
    extprojstack = psi_exit;
    
% if multislice then the exit wave is calculated separately for each tilt angle.  
case 'ms' 
    psi_exit = newim(params2.proc.N, params2.proc.N, nTiltAngles, 'dcomplex');
    %if strcmp(params2.seriesout,'tilt')
    for ll=1:nTiltAngles
        
        tiltang = params2.acquis.tilt(ll);
        fprintf('Simulate tilt angle %3.0f\n', tiltang*180/pi)
        disp(' ')
        disp('Generating Tilt...')
        tic
        [potenext, n] = tiltingMS(poten0,tiltang, params2);
        clear poten0
        toc
        disp(' ')
        sizepot = size(potenext);
        Nm      = max(sizepot(1), sizepot(2));
        thicknessfull = sizepot(3)*voxSz;

        % construct imaginary part 
        if params2.spec.imagpot == 0
            potenfull = potenext;
        elseif params2.spec.imagpot == 1
            potenfull = potenext +1i*params2.spec.potenampl*potenext;
        elseif params2.spec.imagpot == 2 
            potenfull = potenext +1i*params2.spec.potenampl*(newim(potenext)+1);
        elseif params2.spec.imagpot == 3
            imagPot = params2.spec.potenampl*(potenext==0) + params2.spec.proteinampl*(potenext>0);
            potenfull = potenext +1i*imagPot; 
        else
            error('This option for the params.spec.potenampl is not valid. Please choose between 0-3');
        end
        
        clear potenext
        
        
        
        %Fourier domain
        xwm = (voxSz)*Nm;%pixelsize for multislice * size sample
        q_true_pix_m = 1/xwm;
        q_m = rr([Nm Nm])*q_true_pix_m; % frequencies in Fourier domain

        % propagator function Fresnel Propagation 
        dzprop = thicknessfull/n;
        
        % MULTISLICE
        psi_slice = newim(Nm,Nm,'scomplex')+1;

        % Phase grating & Multislice at  same time
        disp('-.-.-.-.-.-.-.-.- STARTING MULTISLICE -.-.-.-.-.-.-.-.-...')
        disp(' ')
        tic
        for ii = 0:n-1
            vai = tic;
                disp('Loading potential...')
                
                
                disp(['SLICE ' char(num2str(ii+1)) '/' char(num2str(n))])
                disp(' ')
                disp('     Phase grating...')
                tic
                % Projected potential in a slice 
                t = mean(potenfull(:,:,ii*(sizepot(3)/n):(ii+1)*(sizepot(3)/n)-1),[],3); 
                toc

                disp(' ')
                disp('     Transmission function...')
                tic
                % transmission functions for a single slice as we go through
                psi_t = exp(1i*params2.inter.sig_transfer*t*dzprop); 
                toc

                disp(' ')
                disp('     Fresnel propagator...')
                tic
                % Fresnel propagator
                P = exp(-1i*pi*params2.inter.lambda*(q_m.^2)*dzprop); % Fresnel propagator
                toc

                disp(' ')
                disp('     Wave propagation...')
                tic
                % Exit wave from current slice
                psi_slice = ift(ft(psi_slice*squeeze(psi_t(:,:,0)))*P);
                toc
            disp(' ')
            disp(['FINISHED SLICE ' char(num2str(ii+1)) '/' char(num2str(n))])
            toc(vai)
        end
        PsiExit = psi_slice;
        psi_exit(:,:,ll-1) = PsiExit;
        disp(' ')
        disp('-.-.-.-.-.-.-.-.- FINISHED MULTISLICE -.-.-.-.-.-.-.-.-')
        toc
        disp(' ')

        
        if strcmp(params2.spec.source, 'amorph')
            thicknessfull = params2.spec.thick/cos(tiltang);
            psi_exit(:,:,ll-1) = psi_exit(:,:,ll-1)*exp(-params2.inter.sig_transfer*params2.spec.potenampl*thicknessfull);
        end
        
        %projpot_r = resample(squeeze(real(psi_exit(:,:,ll-1))),[(voxSz*1e-10/params2.acquis.pixsize) (voxSz*1e-10/params2.acquis.pixsize)]);
        %projpot_i = resample(squeeze(imag(psi_exit(:,:,ll-1))),[(voxSz*1e-10/params2.acquis.pixsize) (voxSz*1e-10/params2.acquis.pixsize)]);
        %proj_pot  = projpot_r+1i*projpot_i; % projected potential
        %%extproj   = extend_exit(squeeze(proj_pot),params2);
        extprojstack(:,:,ll-1) = squeeze(psi_exit(:,:,ll-1));
        
        
    end
end
if strcmp(params2.seriesout,'defocus') || strcmp(params2.seriesout,'dose')
   psi_exit = repmat(psi_exit,[1 1 Nseries]);
end

%% ---------------------------------- CTF with df, ast, envelopes and optionally phase plate

psi_exit=double(psi_exit);   
btot_i=double(btot_i);
for jjj= 1:Nseries
    if strcmp(params2.seriesout,'defocus')
        params2.acquis.df_run=params2.acquis.df(jjj);
    else
        params2.acquis.df_run=params2.acquis.df(1);
    end
    if ~mod(jjj,5)||~mod(jjj,Nseries)
        fprintf(['Calculate the CTF for the ' params2.seriesout sprintf(' series. Image number %3d of %3d\n',  jjj, Nseries)]);
    end
  
  [ctf] = simulateCTF(params2);
  if params2.disp.ctf 
        ctf_out(:,:,jjj)=double(ctf);
  end
  
% For dark field image  
%   central_block = ones(params2.proc.N, params2.proc.N);
%   central_block(((params2.proc.N/2-31):(params2.proc.N/2+32)), ((params2.proc.N/2-31):(params2.proc.N/2+32))) = 0.1 ;  % block the central direct beam 
%   central_block((params2.proc.N/2+1), (params2.proc.N/2+1)) = 0.025 ;
%   
%   Image_spp=dip_image(ones(params2.proc.N,params2.proc.N),'complex');
%         x_ic = xx(Image_com);         
%         y_ic = yy(Image_com);
%         phi_ic = atan2(y_ic,x_ic);   
%         SPhPlate = exp(1i*params2.mic.PP_Phase*phi_ic).*Image_com; % For spiral phase plate
%         SPhPlate = exp(1i*params2.mic.PP_Phase)*Image_spp;
%         SPhPlateIm = gaussf(imag(SPhPlate),1);
%         SPhPlateRe = gaussf(real(SPhPlate),1);
%         SPhPlate = SPhPlateRe+1i*SPhPlateIm;
  
% btot = central_block*SPhPlate*dip_fouriertransform(dip_image(psi_exit(:,:,jjj)),'forward',[1 1 ]); %0]);  % for dark field image, no CTF
  
  
  btot = ctf*dip_fouriertransform(dip_image(psi_exit(:,:,jjj)),'forward',[1 1 ]); %0]);
  btot_i(:,:,jjj) = double(abs(dip_fouriertransform(btot,'inverse',[1 1])).^2); % intensity in the image without camera influence              
end
noiseless_tilt_series = dip_image(btot_i); 


            
%%  --------------------------------- Camera influence          
noiseless_tilt_series=double(noiseless_tilt_series);
series = double(series);
for iii= 1:Nseries
    if strcmp(params2.seriesout,'dose')
        params2.influx=params2.acquis.dose(iii);
    else
        params2.influx=params2.acquis.dose(1);
    end
    if ~mod(iii,5)||~mod(iii,Nseries)
        fprintf(['Calculate the DQE for the ' params2.seriesout sprintf(' series. Image number %3d of %3d\n', iii, Nseries)]);
    end
IntenDetect = DetectSim(squeeze(dip_image(noiseless_tilt_series(:,:,iii))), params2); 
series(:,:,iii) = double(IntenDetect);
end
series = dip_image(series);

 %% ---------------------------------Output structure
imStructOut.series           = series;
imStructOut.noiseless_series = noiseless_tilt_series;
imStructOut.exit             = extprojstack;
if params2.disp.ctf 
   imStructOut.exit = ctf_out;
end
if params2.disp.mtfdqe
   imStructOut.mtf = params2.cam.mtf;
   imStructOut.dqe = params2.cam.dqe;
end