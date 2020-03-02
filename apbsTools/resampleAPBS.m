function PotCub1A = resampleAPBS(filename,wr)
%resampleAPBS: resamples the APBS volume to volume with 1A cubical voxels
% SYNOPSIS:
% PotCub1A = resampleAPBS(filename,wr)
%
% PARAMETERS:
%   filename: apbs output file name
%         wr: flag to wrtie the volume (1 or 0)
%
% OUTPUT:
%   PotCub1A: Output volume with 1A cubical voxels

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

[potential3D, thickness, dxnew, dynew, dznew ] = readinginpotentials3D(filename);
poten = resample(potential3D,[dxnew/dynew dynew/dynew dznew/dynew]); % make the voxels cubic
potential_out1 = resample(poten,[dynew/1 dynew/1 dynew/1]); % resample it to the 1A voxel size
 
nc = phys_const();
Temp = 291; %K =18C % room temperature
convtovolts = nc.kb*Temp/nc.el; % conversion from kT/e to volts (used when calculating PB) 
PotCub1A = potential_out1*convtovolts;

if wr==1
    tom_mrcwrite(single(PotCub1A), 'name', sprintf('%s', [filename(1:end-3), '_apbsCub1A']));  
end
 
 