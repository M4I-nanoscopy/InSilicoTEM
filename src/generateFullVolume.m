function [VolPot, PosOrient] = generateFullVolume(PartPot,params2,circ,mindista)
%generateFullVolume constructs the full specimen volume and places the particles with the volume
% SYNOPSIS:
% [VolPot] = generateFullVolume(PartPot,params2)
%
% PARAMETERS:
%   PartPot: the potentaila of a particle (obsolite for params2.spec.source 'pdb' )
%   params2: structure containing various input paramters 
%      geom: file containing the position (x,y,z relative to the center of the volume) and 
%            Euler orientations of each particle (size: number of particles x 6)
% OUTPUT:
%  VolPot: Interaction potential volume (slab geometry) 

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
szz     = ceil(params2.spec.thick/params2.acquis.pixsize) - mod(ceil(params2.spec.thick/params2.acquis.pixsize),2); % limted by the specified thickness

tic
if strcmp(params2.spec.source, 'pdb')
    
    
    % check the maximal thickness of the particles in case the input is pdb maps 
    mindist  = 0;
    mindistZ = 0;
    dir0 = params2.proc.rawdir;
    if params2.spec.imagpot ~= 3
        list = dir([dir0 filesep 'Particles' filesep params2.spec.pdbin '*MF' sprintf('%3.1f',params2.spec.motblur) '_VoxSize' sprintf('%02.2f',params2.acquis.pixsize*1e10) '*A.raw']);
    else
        list = dir([dir0 filesep 'Particles' filesep params2.spec.pdbin '*MF' sprintf('%3.1f',params2.spec.motblur) '_VoxSize' sprintf('%02.2f',params2.acquis.pixsize*1e10) '*A_Volt' sprintf('%03d',params2.acquis.Voltage/1000) 'kV.raw']);
    end
    numPart = min(size(list,1),numpart); 
    if ~numPart
        error('Cannot find particles to load');
    end
    
    for ss=1:numPart
        fileName = list(ss).name;
        [tokens matchstring] = regexp(fileName,'a(\d+)_Nx(\d+)_Ny(\d+)_Nz(\d+)_Alp(\d+\.?\d*)_Bet(\d+\.?\d*)_Gam(\d+\.?\d*)_MF(\d+\.?\d*)_VoxSize(\d+\.?\d*)A','tokens','match');
        nx = str2num(tokens{1}{2});
        ny = str2num(tokens{1}{3});
        nz = str2num(tokens{1}{4});
        sz12m = max(nx,ny);
        mindist = max(mindist,sz12m);
        mindistZ = max(mindistZ,nz);
    end

    
    if mindistZ>szz
        fprintf(['The thickness is limited by the particle to ' sprintf('%3d nm\n', round(mindistZ*params2.acquis.pixsize*1e9))])
        zrange = 0;
    else   
        zrange = szz/2-mindistZ/2;
    end

    if params2.proc.geom % specific positions
       [xt, yt, zt, alphaNo, betaNo, gammaNo] = PartList(params2); 
        xt=xt'; yt=yt'; zt=zt';
    else
       [xt, yt, zt, numpart] = RndPos(params2,mindista, zrange,circ);
    end  

    transl = [xt, yt, zt];
    szz = max(szz,mindistZ)+1;
%     TransName = sprintf('Transl_%s_xt_yt_zt_Vol%02d_%02d_%02d_NrPart%03d',params2.spec.pdbin,N,N,szz,numPart);
%     save([pwd filesep 'Particles' filesep TransName], 'transl'); 
    volstruct = zeros(N,N,szz);
 
    for ss=1:numPart
        fileName = list(ss).name;
        [tokens matchstring] = regexp(fileName,'a(\d+)_Nx(\d+)_Ny(\d+)_Nz(\d+)_Alp(\d+\.?\d*)_Bet(\d+\.?\d*)_Gam(\d+\.?\d*)_MF(\d+\.?\d*)_VoxSize(\d+\.?\d*)A','tokens','match');
        tt=str2num(tokens{1}{1}); 
        nx=str2num(tokens{1}{2});
        ny=str2num(tokens{1}{3});
        nz=str2num(tokens{1}{4});
        Alp=str2num(tokens{1}{5});
        Bet=str2num(tokens{1}{6});
        Gam=str2num(tokens{1}{7});
        MF=str2num(tokens{1}{8});
        VoxS=str2num(tokens{1}{9});
        alphad(ss)= Alp;
        betad(ss) = Bet;
        gammad(ss)= Gam;
        if  params2.spec.imagpot ~= 3
            outputfilename = sprintf('%s_a%04d_Nx%i_Ny%i_Nz%i_Alp%3.1f_Bet%3.1f_Gam%3.1f_MF%3.1f_VoxSize%02.2fA.raw',params2.spec.pdbin,tt, nx, ny, nz, Alp, Bet, Gam, MF,VoxS);
        else
            outputfilename = sprintf('%s_a%04d_Nx%i_Ny%i_Nz%i_Alp%3.1f_Bet%3.1f_Gam%3.1f_MF%3.1f_VoxSize%02.2fA_Volt%03dkV.raw',params2.spec.pdbin,tt, nx, ny, nz, Alp, Bet, Gam, MF,VoxS, params2.acquis.Voltage/1000);
        end
        OutFileName = [dir0 filesep 'Particles' filesep outputfilename];
        disp(OutFileName);
        fid = fopen(OutFileName, 'r');
        if params2.spec.imagpot~=3
            atompot = fread(fid,nx*ny*nz,'*double');
            fclose(fid);
            atompot4 = reshape(atompot,[nx, ny, nz]);
        else
            atompot = fread(fid,2*nx*ny*nz,'*double');
            fclose(fid);
            atompotF = reshape(atompot,[2*nx, ny, nz]);
            atompotR = atompotF(1:nx, :, :);
            atompotI = atompotF(nx+1:2*nx, :, :);
            atompot4 = atompotR + 1i*atompotI;
        end

        xy0 = [N/2+xt(ss) N/2+yt(ss)]+1;
        xy1 = xy0 + [nx ny] - [1,1];
        
        x   = xy0(1)-floor(nx/2):xy1(1)-floor(nx/2);
        y   = xy0(2)-floor(ny/2):xy1(2)-floor(ny/2);
        z   = round(szz/2-nz/2)+zt(ss) : round(szz/2+nz/2)-1+zt(ss); 
            
        %fprintf('Min Position xyz %5.2f %5.2f %5.2f\n',min(x),min(y),min(z));
        %fprintf('Max Position xyz %5.2f %5.2f %5.2f\n',max(x),max(y),max(z));
        
        volstruct(x,y,z) = atompot4;   
    end
    VolPot = dip_image(permute(volstruct,[2 1 3]),'scomplex');  
    clear volstruct;
elseif strcmp(params2.spec.source, 'map')
    atompot4 = double(PartPot);
    [nx,ny,nz]=size(atompot4);
     mindist  = ceil(max([nx ny nz]));
    if params2.proc.geom % specific positions and orientations
       [xt, yt, zt, alphad, betad, gammad] = PartList(params2); 
       xt=xt'; yt=yt'; zt=zt';
    else
       % random uniform orientation roll-pitch-yaw euler angles
        alpha = 2*pi*rand(1, params2.proc.partNum); 
        beta = acos(1-2*rand(1, params2.proc.partNum));
        gamma = 2*pi*rand(1, params2.proc.partNum);
        alphad = rad2deg(alpha); betad = rad2deg(beta); gammad = rad2deg(gamma);
       [xt, yt, zt, numpart] = RndPos(params2,mindist, 0);
    end 
    szz = max(szz,mindist)+1;
    volstruct = zeros(N,N,szz);

    for ss=1:numpart
        atompot4rot= double(rotation3d(atompot4,'direct', [deg2rad(alphad(ss)) deg2rad(betad(ss)) deg2rad(gammad(ss))] , [], 'bspline', 'zero'));
        nxr = size(atompot4rot,1); nyr = size(atompot4rot,2); nzr = size(atompot4rot,3);
        xy0 = [N/2+xt(ss) N/2+yt(ss)];
        xy1 = xy0 + [nxr nyr] - [1,1];
        x   = xy0(1)-floor(nxr/2):xy1(1) - floor(nxr/2);
        y   = xy0(2)-floor(nyr/2):xy1(2) - floor(nyr/2);
        z   = round(szz/2-nzr/2) + zt(ss) : round(szz/2+nzr/2) - 1 + zt(ss); 
        volstruct(x,y,z) = atompot4rot;
    end
    VolPot = dip_image(permute(volstruct,[2 1 3]),'scomplex');  
    
elseif strcmp(params2.spec.source, 'amorph')
    VolPot =  repmat(PartPot, [1 1 2]);
else
    error('Input not known. Available options for params.spec.source: ''pdb'', ''map'', ''amorph''')
end 

if ~strcmp(params2.spec.source, 'amorph')
   PosOrient = [xt, yt, zt, alphad', betad', gammad'];
else
   PosOrient = 0;
end
toc
 
