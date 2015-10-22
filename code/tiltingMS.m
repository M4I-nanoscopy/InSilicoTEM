function [TiltVol, n] = tiltingMS(InVol,tiltang, params2)   
%tiltingMS Tilts (rotates) the input volume
% SYNOPSIS:
% [TiltVol, n] = tiltingMS(InVol,tiltang, params2) 
%
% PARAMETERS:
%   InVol: input volume 
% titlang: tilt angle
% params2: structure containing various input physical and processing parameters
%
% OUTPUT:
%  TiltVol: Tilted volume
%  n: number of slices for multislice

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic
voxSz = params2.acquis.pixsize;

if tiltang ~= 0
    TiltVol = rotation(InVol, tiltang, 2, 'bspline','zero');
%fix bug :rotation.m function doesn't work with complex input volume. Rotate 2 times (real and imag) or ask Bernd? 
    TiltVol(TiltVol < min(InVol)) = min(InVol);
    TiltVol = cut(TiltVol, [size(InVol,1), size(InVol,2), size(TiltVol,3)]);
else
    TiltVol = InVol;
end

thicknessfull = voxSz*size(InVol,3)/cos(tiltang);
% to ensure the integer number of slices
n          = ceil(thicknessfull/params2.inter.msdz);
if params2.inter.msdz>thicknessfull
    n = 1;
    fprintf('The slice thickness can not be larger than specimen thickness. Use only one slice...\n')
end  

if mod(single(thicknessfull/voxSz),single(n))>0
    intslices       = thicknessfull/voxSz - mod(single(thicknessfull/voxSz),single(n))+ n;
else
    intslices       = thicknessfull/voxSz - mod(single(thicknessfull/voxSz),single(n));
end
%specimenblock   = extend(TiltVol, [size(InVol,1), size(InVol,1), round(thicknessfull/voxSz)],'symmetric','zero',1); 
%TiltVol         = cut(TiltVol, [size(InVol,1), size(InVol,1), intslices]);
TiltVol         = extend(TiltVol, [size(InVol,1), size(InVol,1), intslices],'symmetric',0,1);
%TiltVol         = resample(TiltVol, [1 1 double(intslices/round(thicknessfull/voxSz))]); 
sizepot         = size(TiltVol);
if mod(intslices, sizepot(3))==1
    TiltVol     = extend(TiltVol,     [sizepot(1), sizepot(2), sizepot(3)+1],'symmetric','zero',1);
elseif mod(intslices, sizepot(3))==sizepot(3)
    TiltVol     = cut(TiltVol,    [sizepot(1), sizepot(2), sizepot(3)-1]);
end



     
     

        
