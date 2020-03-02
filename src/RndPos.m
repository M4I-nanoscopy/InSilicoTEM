function [xt, yt, zt, numpart] = RndPos(params2,mindist,zrange)
%RndPos Assigns random translations of the particles within the volume.
% Volume is splitted into the same size regions within which particle
% positon is ranodmized
% 
% SYNOPSIS:
% [xt, yt, zt, numpart] = RndPos(params2, mindist, zrange)
%
% PARAMETERS:
% params2: structure containing various input paramters (e.g. # of particles)
% mindist: miximal size of the particle (minimal distance beetween them)
%  zrange: the paricles are also randomized in z direction depending on their thickness and thickness of the specimen 
%
% OUTPUT:
%       x: translation in x direction [in pixles] from the image center
%       y: translation in y direction [in pixles] from the image center
%       z: translation in z direction [in pixles] from the image center
%
% numpart: number of particles (possibly limited by the field of view)

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

N       = params2.proc.N;
numpart = params2.proc.partNum;
% szz     =ceil(params2.spec.thick/params2.acquis.pixsize)-mod(ceil(params2.spec.thick/params2.acquis.pixsize),2); 
% mindist =0;
% mindistZ=0;
% 
% %list = dir([pwd filesep 'Particles' filesep params2.spec.pdbin '*.raw']);
% list= dir([pwd filesep 'Particles' filesep params2.spec.pdbin '*MF' sprintf('%3.1f',params2.spec.motblur) '_VoxSize' sprintf('%02.2f',params2.spec.voxsize) '*A.raw']);
% numPart = min(size(list,1), numpart);
% for ss=1:numPart
%     fileName = list(ss).name;
%     [tokens matchstring] = regexp(fileName,'a(\d+)_Nx(\d+)_Ny(\d+)_Nz(\d+)_Alp(\d+\.?\d*)_Bet(\d+\.?\d*)_Gam(\d+\.?\d*)_MF(\d+\.?\d*)_VoxSize(\d+\.?\d*)A.raw','tokens','match');
%     nx=str2num(tokens{1}{2});
%     ny=str2num(tokens{1}{3});
%     nz=str2num(tokens{1}{4});
%     sz12m=max(nx,ny);
%     mindist = max(mindist,sz12m);
%     mindistZ = max(mindistZ,nz);
% end
% szz = max(szz, mindistZ);

if numpart> floor((N/(mindist)))^2
    fprintf('The maximum number of particles within this field of view is limited to %d\n', floor((N/(mindist)))^2)
    numpart = floor((N/(mindist)))^2;
end

K     = ceil(sqrt(numpart));
side0 = N/K;
  C=[];
    for yy1=1:K;
        for xx1=1:K;
                coord=(mindist)/2+(side0-mindist)*rand(1,2);
                C=[C; (xx1-1)*side0+coord(1,1), (yy1-1)*side0+coord(1,2)];
        end
    end  
    
xt = floor(C(1:numpart,1)-N/2);
yt = floor(C(1:numpart,2)-N/2);
zt = floor(zrange*rand(numpart,1).*(-1).^(round(rand(numpart,1))));
