% this code reads the potentials from IASA, reads output of APBS (PB potential), 
% resamples PB potential to 1A maps
% and combines it with the potential calcualted via APBS

%Specify the name of IASA and PB input
filenameIASA 
filenamePB 

% ------------ reads IASA potential(should be map with 1A voxel size)
potts_r  = tom_mrcread([filenameIASA]);
potts_r    = dip_image(potts_r.Value);

% ------------ reads PB potential 
PotCub1A = resampleAPBS(filenamePB,1);
potpb_rper = permute(PotCub1A, [1 3 2]);
szpb = size(potpb_rper);
% extends the IASA potentials
pottsext = extend(potts_r, szpb,'symmetric','zero',1); 

% combined potential
potcomb = pottsext + potpb_rper;

% write the potential:
 tom_mrcwrite(single(potcomb), 'name', sprintf('%s', [filenameIASA(1:end-4), '_PB.mrc'])); 




   