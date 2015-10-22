% --------------------------------- InSilicoTEM
% The software simulates TEM image formation of
% (macro)molecules by taking into account the interaction potential, 
% microscope aberrations and detector characteristics.
% Reference: M.Vulovic et al., Journal of Structural Biology, Volume 183, 2013, Pp 19�32
% 
% Milos Vulovic 2013

clc; clear; close all;

% ----------------------- General processing parameters (proc field)
nb = 256; % the size of field viewdisp
params.proc.N               = 2840; %2328; % 2620 the effective bin at tilt angle 80 is 256
params.proc.partNum         = 1;   % Number of particles. 
params.proc.geom            = 1;  % Specify orientation and translation of particles in 'PartList.m' (=1) or generate them randomly (=0) 
params.proc.cores           = 4;  % the numer of workers going to beused for parllel calculation. If you have a very large volume sample (ex.  

% ----------------------- Specimen (spec field)
params.spec.source          = 'pdb';        % Options: 'map', 'pdb', or 'amorph'      
   params.spec.pdbin        = '2GTL_s1' ;      % Required if params.spec.source = 'map' or 'pdb' ( 2GTL, 2WRJ, 1SA0, 1RYP or any other pdb entree)
%    params.spec.mapsample    = '2WRJ_VoxSize1.0A.mrc'; % if input is a map the name should contain the info about voxel size with following convention '_VoxSize'%02.2f'
params.spec.potcontribution = 'iasa';    % Potential type. Options: 'iasa' or 'iasa+pb'(note: for 'iasa+pb' you first need to calcualate the pb potential via apbs. See the manual)
params.spec.motblur         = 0;            % Motion blur in [A]
params.spec.thick           = 60e-09;    %50   % Thickness of the specimen[m]. 
params.spec.imagpot         = 3;            % Amplitude contrast flag. Options: (=0, none) (=1 constant Q) (=2 ice "plasmons") (=3 'plasmons" of ice and protein)         


% ---------------------- Electron-specimen interaction (inter field)
params.inter.type           = 'tpga'; % Options: 'pa' (projection), 'wpoa' (weak-phase), 'pa+wpoa', 'tpga'(thick-phase grating), 'ms'(multislice) 
    params.inter.msdz       = 1e-9; % Approximate thickness of the slice for multislice [m] Required if params.inter.type = 'ms'

       
% ---------------------- Microscope (mic field)
params.mic.Cs          = 2.7e-3;   % Spherical aberration  [m]
params.mic.a_i         = 0.030e-3; % Illumination aperture [rad]
params.mic.C_c         = 2.7e-3;   % Chromatic aberration  [m]
params.mic.deltaE      = 0.7;      % Energy spread of the source [eV]
% ----------------------- aperture
params.mic.diam_obj          = 100e-6;   % Diameter of objective aperture [m]
params.mic.foc               = 2.89e-3;   % Focal distance [m]
% ----------------------- ideal phase plate (optional)
params.mic.PPflag           = 0;        % Phase plate flag 
params.mic.PP_Phase         = 0.5*pi;   % Phase plate phase shift [rad]
params.mic.PP_qcuton        = 1/250*1e10;% Cut-on frequency [1/m]. Typical values: 1/25, 1/75, 1/125, 1/250, 1/250


%----------------------- Acquisition settings (acquis field)
params.acquis.pixsize        = 0.25e-9;   % Pixel size in the specimen plane [m]3748
params.acquis.df             = [5000]*1e-9; % Defocus [m]. Note: undefocus df>0; overfocus df <0. 
params.acquis.ast            = [0]*1e-9;    % Astigmatism [m]. 
params.acquis.astangle       = 0*pi/180;    % Astigmatism angle [rad]
params.acquis.Voltage        = 200e3;       % Acceleration voltage
params.acquis.tilt           = [-80:20:0]/180*pi;  % Tilt geometry  
params.acquis.dose_on_sample = 1; %[40]/length(params.acquis.tilt); % Integrated flux (for the full tilt series) [e-/A2]



% ----------------------- Detector-Camera (cam field)
params.cam.type              = 'FalconII'; % Options: 'custom', 'Eagle4k', 'US4000', 'US1000GIF', 'FalconI', 'perfect' (64% at Nq -counting mode), 'ideal' (100% at Nq)
params.cam.bin               = 1; % hardware binning
% if params.cam.type = 'custom' please characterize your camera (provide MTF and DQE files and if needed readnoise and dark current).
% The characterization can be done via e.g. the tools and methods described in Vulovic et al. 2010 Acta Crist. D 
% If you want to simulate MTF without accurate characterization set params.cam.GenMTFasEMG
params.cam.GenMTFasEMG      = 1;
params.cam.DQEflag          = 1; % Flag for dqe (=0 means NTF = MTF)


% ---------------------- Display what? (disp field)
params.disp.generateWhat   = 'im'; % Options: 'im', 'exitw', 'imNoiseless'
params.disp.ctf            = 0; % Flag to display CTF
params.disp.mtfdqe         = 1; % Flag to display MTF and DQE
%----------------------- Delei: parameter
params.particle.rotate = 0; % rotage particle before packing in the full model


% ---------------------- Parse parameters
params2 = parsePar(params);

tic;
%% ---------------------- Generate and/or load 3D potential of a particle
[PartPot, params2] = loadSamples(params2);

fprintf('*******  Sample Load: %4.2f min\n', toc/60);
%% ---------------------- Padding and placing the particles within the volume
tic
[InputVol, PosOrient] = generateFullVolume(PartPot,params2);
fprintf('*******  Generate Volume: %4.2f min\n', toc/60); 

InputVol = InputVol(:,:,1:end); % 

%% {
% ---------------------- Image simulation
tic
[imStructOut] = simTEM(InputVol, params2);
fprintf('******* simTEM: %4.2f min\n', toc/60); 

%%  crop the image size to fit the size of FOV
nb = nb;
proj.exit = cut(imStructOut.exit,nb);
proj.blur = cut(imStructOut.noiseless_series,nb);
proj.finl = cut(imStructOut.series,nb);


% ---------------------- Display
switch params.disp.generateWhat
    case 'im'
%            dipshow(imStructOut.series, 'percentile');
           dipshow(proj.finl, 'percentile');
    case 'exitw'
%            dipshow(imStructOut.exit, 'lin');
           dipshow(proj.exit, 'lin');
    case 'imNoiseless'
%            dipshow(imStructOut.noiseless_series, 'lin');  
           dipshow(proj.blur, 'lin'); 
end

%}
% --------------------- Writing the image (stack in case of a series).      

%     tom_mrcwrite(double(imStructOut.series), 'name','FinalImage')
%     tom_mrcwrite(double(imStructOut.exit), 'name','ExitWave')
%     tom_mrcwrite(double(imStructOut.noiseless_series), 'name','Noiseless')