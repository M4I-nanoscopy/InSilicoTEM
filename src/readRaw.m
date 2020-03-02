% An example of reading raw files.
% this is incorporated in the function "generateFullVolume"

pdbin   = '1SA0';
dirPart = [pwd filesep 'Particles']; 
        fileName = '1SA0_a0001_Nx40_Ny45_Nz24_Alp222.4_Bet127.7_Gam65.9_MF0.0_VoxSize4.00A.raw';
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
        outputfilename = sprintf('%s_a%04d_Nx%i_Ny%i_Nz%i_Alp%3.1f_Bet%3.1f_Gam%3.1f_MF%3.1f_VoxSize%02.2fA.raw',pdbin,tt, nx, ny, nz, Alp, Bet, Gam, MF,VoxS);
        OutFileName = [dirPart filesep outputfilename];
        disp(OutFileName);
        fid = fopen(OutFileName, 'r');
        atompot = fread(fid,nx*ny*nz,'*double');
        fclose(fid);
        atompot4 = reshape(atompot,[nx, ny, nz]); 
        dipshow(dip_image(atompot4))
     