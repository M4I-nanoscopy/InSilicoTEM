function [PartPot, pout] = loadSamples(params2)
%loadSamples Loads or calculates the interaction potential of the particles
%depending on the required input params2.spec.source : 'pdb', 'map' or 'amorph'
% SYNOPSIS:
% [PartPot,pout] = loadSamples(params2)
%
% PARAMETERS:
%  params2: structure containing various input paramters
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

       
pout = params2;

if strcmp(params2.spec.source, 'pdb')   % generate the potential maps from pdb
    % first check if the particles are already generated
    if params2.spec.imagpot~=3
        list = dir([pwd filesep 'Particles' filesep params2.spec.pdbin '*MF' sprintf('%3.1f',params2.spec.motblur) '_VoxSize' sprintf('%02.2f',params2.acquis.pixsize*1e10) '*A.raw']);
    else
        list = dir([pwd filesep 'Particles' filesep params2.spec.pdbin '*MF' sprintf('%3.1f',params2.spec.motblur) '_VoxSize' sprintf('%02.2f',params2.acquis.pixsize*1e10) '*A_Volt' sprintf('%03d',params2.acquis.Voltage/1000) 'kV.raw']);
    end
   params2.NumGenPart = size(list,1);
   if params2.proc.partNum <= params2.NumGenPart %&& ~params2.proc.geom 
         PartPot = 0;
         fprintf('No need for generation of new particles. Load the potential maps from the folder ''Particles''...\n')
   else         
        if params2.proc.geom 
         % load the specified orientation of the particles
         [xt, yt, zt, alphad, betad, gammad] = PartList(params2);
        else
         % random uniform orientation roll-pitch-yaw euler angles
        alpha = 2*pi*rand(1, params2.proc.partNum); 
        beta  = acos(1-2*rand(1, params2.proc.partNum));
        gamma = 2*pi*rand(1, params2.proc.partNum);
        alphad = rad2deg(alpha); betad = rad2deg(beta); gammad = rad2deg(gamma);
        end  
        % calculate the atomic potential from the pdb with given
        % orientation. Save the particles into subfolder 'Particles' (wr=1)
      wr = 1;
      [atompotF, numpart] = AtomPotRot(params2, alphad, betad, gammad, wr);
      pout.proc.partNum = numpart + params2.NumGenPart;
      PartPot = atompotF;
   end
    
elseif strcmp(params2.spec.source, 'amorph') % the specimen is amorphous
    projected = dip_image(randn(params2.proc.N,params2.proc.N));
    DMft      = ft(projected);
    if strcmp(params2.spec.imagpot2specm, 'amorIce')
        d0 = 1.54;% Ang
    elseif strcmp(params2.spec.imagpot2specm, 'amorC')
        d0 = 2.88; % Ang
    end
    dmin = 2*d0/(params2.acquis.pixsize/1e-10);
    Btot = 2*pi^2*(dmin)^2; 
    PartPot = real(ift(DMft*exp(-Btot*(rr(DMft,'freq')/params2.acquis.pixsize*1e-10).^2)));
    
elseif strcmp(params2.spec.source, 'map')  % load already existing maps from folder 'pot'    
    Potdir = [pwd filesep 'MAPs'];
    if ~exist(Potdir, 'dir')
        error('Could not find a suitable folder for potential files. Make a subfolder pot in the working directory or specify the directory in loadSamples.m'); %change from error to message 
    end
    
    if exist([Potdir filesep params2.spec.mapsample]) 
         % load already existing map
         atompot = tom_mrcread([Potdir filesep params2.spec.mapsample]);
         atompot = dip_image(atompot.Value);
         solvpot = mean(atompot(end-2:end,:,end-2:end)); % mean value of the solvent
         atompot = atompot-solvpot;
         [tokens matchstring] = regexp(params2.spec.mapsample,'_VoxSize(\d+\.?\d*)A','tokens','match');
         pout.spec.voxsize = str2num(tokens{1}{1});
         % it is neccesary to apply low-pass filter before eventual downsampling
        if pout.spec.voxsize*1e-10 <= params2.acquis.pixsize
           atompotBl  = gaussf(mat2im(atompot), sqrt((params2.acquis.pixsize/(pout.spec.voxsize*1e-10))^2-1), 'best');
        else
           warning('The voxel size of the available potential is larger than the pixel size')
        end
         PartPot = resample(atompotBl, pout.spec.voxsize*1e-10/params2.acquis.pixsize); 
         if params2.spec.motblur~=0
            PartPot = motionBlur(PartPot,params2);
         end
    else 
         error('The map or its path is not available')
    end     
else
    error('This type of the input type is not known. Options are ''pdb'' ''map'' or '' amorphous''')
end  



