% --------------------------------- TEMsim
% Function based on the original InSilicoTEM to run the simulator from
% external script.
%
% The software simulates TEM image formation of
% (macro)molecules by taking into account the interaction potential, 
% microscope aberrations and detector characteristics.
% Reference: M.Vulovic et al., Journal of Structural Biology, Volume 183, 2013, Pp 19�32
% 
% Milos Vulovic 2013
clc; clear; close all;

% ----------------------- Parameters needed to generate the particles coordinates (proc field)
params.proc.N               = 512;         % Image size (field of view)
params.proc.geom            = 0;            % Specify orientation and translation of particles in 'PartList.m' (=1) or generate them randomly (=0) 
params.proc.cores           = 24;           % the numer of matlab pools to be open for parfor loops 
params.proc.rawdir          = 'Raw';      % Specify the directory to save the particles' potential files.

% ----------------------- Specimen (spec field)
params.spec.source           = 'pdb';                   % Options: 'map', 'pdb', or 'amorph'      
params.spec.pdbin            = ["2fha"];  % Put the model which needs to be at the center of the micrograph at the end.
                                                        % Required if params.spec.source = 'map' or 'pdb' ( 2GTL, 2WRJ, 1SA0, 1RYP or any other pdb entree)
params.spec.mapsample        = '2WRJ_VoxSize1.0A.mrc';  % if input is a map the name should contain the info about voxel size with following convention '_VoxSize'%02.2f'
params.spec.potcontribution  = 'iasa';                  % Potential type. Options: 'iasa' or 'iasa+pb'(note: for 'iasa+pb' you first need to calcualate the pb potential via apbs. See the manual)
params.spec.motblur_series   = [0];                   % Motion blur in [A], (If multiple MB, enter them as a vector)
params.spec.imagpot          = 3;                       % Amplitude contrast flag. Options: (=0, none) (=1 constant Q) (=2 ice "plasmons") (=3 'plasmons" of ice and protein)         
params.spec.partdistribution = '3D';                    % 2D or 3D Fill the volumn with particles of 2D nonoverlapped or 3D nonoverlapped.
params.spec.thick            = [200]*1e-10;             % (Fill the blank in A) Thickness of the specimen[m]. 
params.spec.partdensity      = 1;                       % (0,1] Density of the particles, percentage of the maximum number.

% ---------------------- Electron-specimen interaction (inter field)
params.inter.type           = 'wpoa';       % Options: 'pa' (projection), 'wpoa' (weak-phase), 'pa+wpoa', 'tpga'(thick-phase grating), 'ms'(multislice) 
params.inter.msdz           = 2e-9;       % Approximate thickness of the slice for multislice [m] Required if params.inter.type = 'ms'

       
% ---------------------- Microscope (mic field)
params.mic.Cs          = 2.7e-3;            % Spherical aberration  [m]
params.mic.a_i         = 0.030e-3;          % Illumination aperture [rad]
params.mic.C_c         = 2.7e-3;            % Chromatic aberration  [m]
params.mic.deltaE      = 0.7;               % Energy spread of the source [eV]

% ----------------------- aperture
params.mic.diam_obj          = 100e-6;      % Diameter of objective aperture [m]
params.mic.foc               = 4.7e-3;      % Focal distance [m]

% ----------------------- ideal phase plate (optional)
params.mic.PPflag             = 0;         % Phase plate flag 
    params.mic.PP_Phase       = pi/2;       % Phase plate phase shift [rad] *pi
    params.mic.PP_qcuton      = 1/128*1e10; % Cut-on frequency [1/m]. Typical values: 1/25, 1/75, 1/125, 1/250, 1/250

params.mic.SPPflag            = 0 ;         % Spiral Phase Plate Flag
    params.mic.SPP_Phase      = 1;          % quantum number l
    params.mic.SPP_qcuton     = 1/128*1e10;  % Cut-on frequency [1/m]. Typical values: 1/25, 1/75, 1/125, 1/250, 1/250

%----------------------- Acquisition settings (acquis field)
defocus_range                       = [1200 1200];    % Defocus range [nm]
                                                     % Blank line reversed for defocus
params.acquis.pixsize               = [1]*1e-10;    % (Fill the blank in A) Pixel size in the specimen plane [m]3748
params.acquis.ast                   = 0*1e-9;        % Astsaved = sigmatism [m]. 
params.acquis.astangle              = 0*pi/180;      % Astigmatism angle [rad]
params.acquis.Voltage               = 300e3;         % Acceleration voltage
params.acquis.tilt                  = [45]/180*pi;      % Tilt geometry  
params.acquis.dose_on_sample_series = [40]/length(params.acquis.tilt); % Integrated flux (for the full tilt series) [e-/A2], if multiple flux, enter as a vector

% ----------------------- Detector-Camera (cam field)
params.cam.type              = 'FalconIII_EC'; % Options: 'custom', 'Eagle4k', 'US4000', 'US1000GIF', 'FalconI', 'FalconIII_Linear','FalconIII_EC', 'perfect' (64% at Nq -counting mode), 'ideal' (100% at Nq)
params.cam.bin               = 1; % hardware binning
% if params.cam.type = 'custom' please characterize your camera (provide MTF and DQE files and if needed readnoise and dark current).
% The characterization can be done via e.g. the tools and methods described in Vulovic et al. 2010 Acta Crist. D 
% If you want to simulate MTF without accurate characterization set params.cam.GenMTFasEMG
params.cam.GenMTFasEMG      = 1;
params.cam.DQEflag          = 1; % Flag for dqe (=0 means NTF = MTF)
params.cam.cf_series        = [14]; % Conversion factor,(If multiple CF, enter them as a vector).
                                    % This is for FalconIII at 300 keV, need to specify for other camera.
                                    % Modified with Rado 1.62 A data from EMPAIR.                               

% ---------------------- Display what? (disp field)
params.disp.generateWhat   = 'im'; % Options: 'im', 'exitw', 'imNoiseless'
params.disp.ctf            = 1; % Flag to display CTF
params.disp.mtfdqe         = 0; % Flag to display MTF and DQE

% Run simulator under different conditions:                               %
%       - defocus range                                                   %
%       - motion blur (radiation damage)                                  %
%       - correction factor of the microscope                             %
%       - electron dose                                                   %
%       - size of the micrograph (pixel number)                           %
%       - pixel size of the micrograph                                    %
%       - distance between particles                                      %

% More parameters to be defined:
dir0 = '/mnt/turbo/y.zhang';  % Select folder where to save micrographs
workstation = 'giga';                  % Mark workstation (avoid name conflict)
mg = 1;                                 % Number of micrographs to be generated
mindist = (160*1e-10)/(params.acquis.pixsize);    % MiniImum distance between particles divided by pixsize, 
                                                          % Depends on type of protein (apo ferrtin ~150)
count = 0;
tot = mg*length(params.spec.motblur_series)*length(params.cam.cf_series)*length(params.acquis.dose_on_sample_series); % Total number of micrographs
rng('shuffle');

time = tic;
for micro = 1: mg
    for cf = 1 : length(params.cam.cf_series)
        params.cam.cf = params.cam.cf_series(cf);
        for mb = 1 : length (params.spec.motblur_series)
            params.spec.motblur = params.spec.motblur_series(mb);
            for d = 1:length(params.acquis.dose_on_sample_series)
                params.acquis.dose_on_sample = params.acquis.dose_on_sample_series(d);
                
                if isequal(params.acquis.dose_on_sample, params.acquis.dose_on_sample_series(1))    
                % for focal pair, particles at same position on micrographs
                % at differentfloor(p/6);floor(p/6);floor(p/6);floor(p/6);floor(p/6) defocus value
%                 try
                count = count +1;
                disp(...
                [char("###################### Starting new Simulation ######################") newline...
                 char("       Micrographs: " + tot) newline...
                 char("       Size:        " + params.proc.N + "x" + params.proc.N ) newline...
                 char("       Particles:   " + "Randomized") newline...
                 char("       Defocus:     [" + defocus_range(1) + " - " + defocus_range(2)+"]nm") newline...
                 char("       Motion blur: " + params.spec.motblur(1)) newline...
                 char("       Corr factor: " + params.cam.cf(1)) newline...
                 char("       Dose:        " + params.acquis.dose_on_sample(1)) 'e/A^2' newline...
                 '                                 ' char(datetime(now,'ConvertFrom','datenum'))])             
                
                s0 = 'Micrographs%04d';
                s01 = sprintf(s0, micro);
                partspec = sprintf('%s_%s.txt',s01, workstation);
                params.spec.mgraphsNum      = partspec; % file records the particles translational and rotational parameters                
                
                % Randomize particles positions for each micrograph
                [circl, p] = GetParticleCoordinate(mindist, params);
                params.proc.partNum = [p]; % Number of particles, the sequence of number should match the particle's type respectively.               
                
                % Randomize defocus for each micrographs
                defocus = (defocus_range(2)-defocus_range(1))*rand(1,1)+defocus_range(1);
                params.acquis.df = [defocus]*1e-9;       % Defocus range[m]. Note: undefocus df>0; overfocus df <0. 
                
                save('positionAnddefocus.mat', 'circl', 'defocus');
                
                disp([newline newline newline])                  
                disp( "-----------------------   Progress:   "+count+" / "+tot+" Micrographs   -----------------------")
                disp(' ')
                disp(char( "Particles:    "+params.proc.partNum) )
                disp([char( "Defocus:      "+defocus) 'nm' newline])
                disp(' ')
                
                % ---------------------- Parse parameters
                params2 = parsePar(params);
                % ---------------------- Generate and/or load 3D potential of a particle
                [PartPot, params2] = loadSamples(params2);
                % ---------------------- Padding and placing the particles within the volume
                if strcmp(params.inter.type,'wpoa')
                    [InputVol] = generateFullVolume_wpoa(PartPot,params2,circl,mindist);
                else
                    [InputVol] = generateFullVolume(PartPot,params2,circl,mindist);
                end
                % ---------------------- Image Formation---------------------------------------------------------------------------------------------------------           
                [imStructOut] = simTEM(InputVol, params2);
                % ---------------------- Display
%                 switch params.disp.generateWhat
%                     case 'im'
%                         dipshow(imStructOut.series, 'lin')
%                     case 'exitw'
%                         dipshow(imStructOut.exit, 'lin')
%                     case 'imNoiseless'
%                         dipshow(imStructOut.noiseless_series, 'lin')           
%                 end

                % ---------------------- Display noiseless
                % for i = 1:size(imStructOut.noiseless_series,3)   
                %     if size(imStructOut.noiseless_series,3) > 1
                %         subplot(3,3,i)
                %     else
                %         %figure
                %     end
                % %     imagesc(imStructOut.noiseless_series(:,:,i))
                %     dipshow(imStructOut.noiseless_series, 'lin')  
                % end
                % colormap gray
                ImageOut = imStructOut.series;
%                delete ./Raw/Particles/*.raw;  % delete particle positions for the generated micrographs (for memory)

                % Write micrograph to .MRC file
                s1 = "MicrographNr"+ (micro) +"_size"+params.proc.N+"x"+params.proc.N+"_pixsize"+params.acquis.pixsize*1e12+ ...
                    "_partnr"+p+"_dose"+params.acquis.dose_on_sample_series(d)+"_cf"+params.cam.cf_series(cf)+"_mb"+ ...
                    params.spec.motblur_series(mb)+"_df"+round(defocus)+"_PP"+params.mic.PPflag+".mrc";
                s11 = sprintf('%s_%s', workstation, s1);
                if ~ exist([dir0 filesep 'Micrographs'])
                    mkdir([dir0 filesep], 'Micrographs');
                end
                dir1 = '/mnt/turbo/y.zhang/Micrographs';
                WriteMRC(double(ImageOut.'),1.0,[dir1 filesep char(s11)]);
                WriteMRC(double(imag(imStructOut.exit)),1.0,[dir1 filesep char("CTF_"+char(s11))]);
                disp(' ')
                disp('Successful Micrograph')
                
                % Write tiltseries to .avi file
                if length(params.acquis.tilt) > 1
                    name = sprintf('%s.avi',s11);
                    if ~ exist([dir0 filesep 'Video'])
                        mkdir([dir0 filesep], 'Video');
                    end
                    dir2 = '/mnt/ssd/y.zhang/Video';
                    video = VideoWriter(sprintf('%s/%s',dir2,name));      %create the video object
                    open(video);                    %open the file for write

%                     for ss = 1: length(params.acquis.tilt)
%                         I = mat2gray(double(ImageOut(:,:,ss-1)));
%                         writeVideo(video,I);        % write the image to file
%                     end
%                     close(video);
                end
                    
%                 catch ex
%                     % In case the simulator fails to produce a micrograph,
%                     % it will skip it and display a warning
%                       warning(['Failed Micrograph' newline ex.identifier  ...
%                                     newline ex.message      ...
%                                     newline ex.stack.name])
%                       disp(' ')
%                       
%                 end
                
                end
            end
        end
    end
end

disp(' ')
disp(...
   [char("###################### End of Simulation ######################") newline newline...
    '                     ' char(datetime(now,'ConvertFrom','datenum')) newline])
toc(time)
                
