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
N = params2.proc.N;
InVol = dip_array(InVol);
if tiltang ~= 0
    sizeVol = size(InVol);
    TiltVol = zeros(floor(sizeVol(1)*cos(abs(tiltang))+sizeVol(3)*sin(abs(tiltang))),sizeVol(2), ...
        floor(sizeVol(1)*sin(abs(tiltang))+sizeVol(3)*cos(abs(tiltang))), 'single');
    for ii = 1:sizeVol(1)
            for kk = 1:sizeVol(3)
                if tiltang > 0
                    TiltVol(round(ii*cos(tiltang)+kk*sin(tiltang)),:,...
                        round(-ii*sin(tiltang)+kk*cos(tiltang)+abs(-sizeVol(1)*sin(tiltang)+cos(tiltang)))+1) = InVol(ii,:,kk);
                else 
                    TiltVol(round(ii*cos(tiltang)+kk*sin(tiltang)+abs(cos(tiltang)+sizeVol(3)*sin(tiltang)))+1,:,...
                        round(-ii*sin(tiltang)+kk*cos(tiltang))) = InVol(ii,:,kk);
                end
            end
    end
    clear InVol;
    TiltVol = TiltVol(floor((size(TiltVol,1)+2-N)/2):floor((size(TiltVol,1)+N)/2),floor((size(TiltVol,2)+2-N)/2):floor((size(TiltVol,2)+N)/2),:);
%     TiltVol = dip_image(TiltVol);
else  
    TiltVol = InVol;
    clear InVol;
    TiltVol = TiltVol(floor((size(TiltVol,1)+2-N)/2):floor((size(TiltVol,1)+N)/2),floor((size(TiltVol,2)+2-N)/2):floor((size(TiltVol,2)+N)/2),:);
end

thicknessfull = voxSz*size(TiltVol,3);
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
TiltVol         = padarray(TiltVol, [0, 0, double(intslices-thicknessfull/voxSz)],'post');
TiltVol = dip_image(TiltVol);
%TiltVol         = resample(TiltVol, [1 1 double(intslices/round(thicknessfull/voxSz))]); 
sizepot         = size(TiltVol);
if mod(intslices, sizepot(3))==1
    TiltVol     = extend(TiltVol,     [sizepot(1), sizepot(2), sizepot(3)+1],'symmetric','zero',1);
elseif mod(intslices, sizepot(3))==sizepot(3)
    TiltVol     = cut(TiltVol,    [sizepot(1), sizepot(2), sizepot(3)-1]);
end