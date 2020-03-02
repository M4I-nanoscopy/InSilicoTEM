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

function [imStructOut] = TEMsim(df,mb,cf,saved,rpdb,p,dose,circl,pix,pixsize,pp,mindist)
% ----------------------- General processing parameters (proc field)
params.proc.N               = pix;      % Image size (field of view)
params.proc.partNum         = p;            % Number of particles. 
params.proc.geom            = 0;            % Specify orientation and translation of particles in 'PartList.m' (=1) or generate them randomly (=0) 
params.proc.cores           = 10;           % the numer of matlab pools to be open for parfor loops 
params.proc.rawdir          = './Raw';


% ----------------------- Specimen (spec field)
params.spec.source          = 'pdb';                % Options: 'map', 'pdb', or 'amorph'      
params.spec.pdbin        = '1dpx' ;               % Required if params.spec.source = 'map' or 'pdb' ( 2GTL, 2WRJ, 1SA0, 1RYP or any other pdb entree)
params.spec.mapsample    = '2WRJ_VoxSize1.0A.mrc';  % if input is a map the name should contain the info about voxel size with following convention '_VoxSize'%02.2f'
params.spec.potcontribution = 'iasa';               % Potential type. Options: 'iasa' or 'iasa+pb'(note: for 'iasa+pb' you first need to calcualate the pb potential via apbs. See the manual)
params.spec.motblur         = mb;                   % Motion blur in [A]
params.spec.thick           = 10e-09;               % Thickness of the specimen[m]. 
params.spec.imagpot         = 0;                    % Amplitude contrast flag. Options: (=0, none) (=1 constant Q) (=2 ice "plasmons") (=3 'plasmons" of ice and protein)         


% ---------------------- Electron-specimen interaction (inter field)
params.inter.type           = 'ms';       % Options: 'pa' (projection), 'wpoa' (weak-phase), 'pa+wpoa', 'tpga'(thick-phase grating), 'ms'(multislice) 
    params.inter.msdz       = 2e-9;         % Approximate thickness of the slice for multislice [m] Required if params.inter.type = 'ms'

       
% ---------------------- Microscope (mic field)
params.mic.Cs          = 2.7e-3;            % Spherical aberration  [m]
params.mic.a_i         = 0.030e-3;          % Illumination aperture [rad]
params.mic.C_c         = 2.7e-3;            % Chromatic aberration  [m]
params.mic.deltaE      = 0.7;               % Energy spread of the source [eV]
% ----------------------- aperture
params.mic.diam_obj          = 100e-6;      % Diameter of objective aperture [m]
params.mic.foc               = 4.7e-3;      % Focal distance [m]
% ----------------------- ideal phase plate (optional)
params.mic.PPflag             = pp;         % Phase plate flag 
    params.mic.PP_Phase       = pi/2;       % Phase plate phase shift [rad] *pi
    params.mic.PP_qcuton      = 1/128*1e10; % Cut-on frequency [1/m]. Typical values: 1/25, 1/75, 1/125, 1/250, 1/250

params.mic.SPPflag            = 0 ;         % Spiral Phase Plate Flag
    params.mic.SPP_Phase      = 1;          % quantum number l
    params.mic.SPP_qcuton     = 1/128*1e10;  % Cut-on frequency [1/m]. Typical values: 1/25, 1/75, 1/125, 1/250, 1/250

%----------------------- Acquisition settings (acquis field)
params.acquis.pixsize        = (pixsize)*1e-10;       % Pixel size in the specimen plane [m]3748
params.acquis.df             = df*1e-9;       % Defocus [m]. Note: undefocus df>0; overfocus df <0. 
params.acquis.ast            = 0*1e-9;        % Astsaved = 0igmatism [m]. 
params.acquis.astangle       = 0*pi/180;      % Astigmatism angle [rad]
params.acquis.Voltage        = 300e3;         % Acceleration voltage
params.acquis.tilt           = 0/180*pi;      % Tilt geometry  
params.acquis.dose_on_sample = [dose]/length(params.acquis.tilt); % Integrated flux (for the full tilt series) [e-/A2]

% ----------------------- Detector-Camera (cam field)
params.cam.type              = 'FalconIII_EC'; % Options: 'custom', 'Eagle4k', 'US4000', 'US1000GIF', 'FalconI', 'FalconIII_Linear','FalconIII_EC', 'perfect' (64% at Nq -counting mode), 'ideal' (100% at Nq)
params.cam.bin               = 1; % hardware binning
% if params.cam.type = 'custom' please characterize your camera (provide MTF and DQE files and if needed readnoise and dark current).
% The characterization can be done via e.g. the tools and methods described in Vulovic et al. 2010 Acta Crist. D 
% If you want to simulate MTF without accurate characterization set params.cam.GenMTFasEMG
params.cam.GenMTFasEMG      = 1;
params.cam.DQEflag          = 1; % Flag for dqe (=0 means NTF = MTF)


% ---------------------- Display what? (disp field)
params.disp.generateWhat   = 'im'; % Options: 'im', 'exitw', 'imNoiseless'
params.disp.ctf            = 0; % Flag to display CTF
params.disp.mtfdqe         = 0; % Flag to display MTF and DQE



% ---------------------- Parse parameters
params2 = parsePar(params,cf);


% ---------------------- Generate and/or load 3D potential of a particle

if saved
    a = rpdb;
else
    a = 0;
end

[PartPot, params2] = loadSamples(params2,a,saved);


% ---------------------- Padding and placing the particles within the volume
[InputVol, PosOrient] = generateFullVolume(PartPot,params2,circl,mindist);


% ---------------------- Image Formation
[imStructOut] = simTEM(InputVol, params2);


% ---------------------- Display
switch params.disp.generateWhat
    case 'im'
           dipshow(imStructOut.series, 'lin')
    case 'exitw'
           dipshow(imStructOut.exit, 'lin')
    case 'imNoiseless'
           dipshow(imStructOut.noiseless_series, 'lin')           
end

% ---------------------- Display noiseless
% for i = 1:size(imStructOut.noiseless_series,3)
%     
%     
%     if size(imStructOut.noiseless_series,3) > 1
%         subplot(3,3,i)
%     else
%         %figure
%     end
% %     imagesc(imStructOut.noiseless_series(:,:,i))
% 
%     dipshow(imStructOut.noiseless_series, 'lin')  
% 
% end
% colormap gray



end