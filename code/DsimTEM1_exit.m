function psi_exit = DsimTEM1_exit(InputVol,params2)

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
extprojstack  = newim(params2.proc.N,params2.proc.N,Nseries, 'complex');

% btot_i        = newim(params2.proc.N,params2.proc.N,Nseries);
% series        = newim(params2.proc.N,params2.proc.N,Nseries);
poten0 = InputVol;
thickness = voxSz*size(poten0,3);

switch params2.inter.type
    case{'pa+wpoa', 'pa', 'wpoa', 'tpga'}
        
        poten0 = permute(poten0,[1 3 2]);
        % calculate the projected potential (neglects the thickness of the specimen)
        switch params2.inter.type
            case {'pa+wpoa','pa'}
                poten0_stack = FreqDomRadon_2d(poten0,params2.acquis.tilt);
                % in case of wpoa and tpga the Ewald sphere (parabola) is sampled and propagation through the specimen thickness is taken into account
            case {'wpoa','tpga'}
                poten0_stack = FreqDomRadon_Ewald_2d(poten0,params2.acquis.tilt, params2);
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
            [potenext, n] = tiltingMS(poten0,tiltang, params2);
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
            
            % preallocation for multislice
            t  = newim([Nm, Nm, n],'dcomplex');
            
            % projected potential within slices (phase grating)
            for ii = 0:n-1
                t(:,:,ii) = mean(potenfull(:,:,ii*(sizepot(3)/n):(ii+1)*(sizepot(3)/n)-1),[],3);
            end
            
            xwm = (voxSz)*Nm;%pixelsize for multislice * size sample
            %Fourier domain
            q_true_pix_m = 1/xwm;
            q_m = rr([Nm Nm])*q_true_pix_m; % frequencies in Fourier domain
            
            % propagator function Fresnel Propagation
            dzprop = thicknessfull/n;
            % transmission functions; the slice thickness is constant
            psi_t = exp(1i*params2.inter.sig_transfer*t*dzprop);
            % multislice
            PsiExit = multislice(psi_t, Nm, n, params2.inter.lambda, q_m, dzprop);
            psi_exit(:,:,ll-1) = PsiExit;
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