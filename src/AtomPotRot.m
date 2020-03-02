function [outpot,numpartExtra,rpdb] = AtomPotRot(params2, alphad, betad, gammad, wr,rpdbsave,saved,process)
%AtomPotRot Calculates the interaction potential of the particles. Inspired
% by the pdb.c of TEMSimulator (H. Rullg�rd et al, Journal of Microscopy 243 (3) (2011)) 
% SYNOPSIS:
% [outpot,numpartExtra] = AtomPotRot(params2, alphad, betad, gammad, wr)
%
% PARAMETERS:
%   params2: structure containing various input paramters (e.g. input pdb file, voxel size)
%   alphad : the first Euler angles (in degrees)
%    betad : the second Euler angles (in degrees)
%   gammad : the third Euler angles (in degrees)
%       wr : writing out the potential map in folder 'Particles'
%
% OUTPUT:
%        outpot: potential map of only the last particle (obsolite). All the potentials are saved in 'Particles' folder
%  numpartExtra: number of particle generated



pdbin      = params2.spec.pdbin;
voxel_size = params2.spec.voxsize;
Vwat       = params2.spec.Vwat;
pixsz      = params2.acquis.pixsize*1e10;
alpha = deg2rad(alphad); beta = deg2rad(betad); gamma = deg2rad(gammad);
dir0 = params2.proc.rawdir;

% Read the PDB file
wrPdb = sprintf('%s.pdb',pdbin);
%wrMrc = sprintf('%s.mrc',pdbin);
path2pdb = [pwd filesep 'PDBs' filesep wrPdb];
if ~exist(path2pdb)
    getpdb(pdbin,'TOFILE',path2pdb);
end

if ~saved
    tic
    disp('Loading PDB file...')
    rpdb = pdbread(path2pdb);
    toc
else
    disp('Using pre-saved PDB file...')
    rpdb = rpdbsave;
end

assignin('base','rpdbsave',rpdb);


oversam  = 1; % oversampling of the potential 4 means map with voxel size of 0.25 A.
extrapix = 10; % extra voxels to extend the volume 10A = 1nm

% extract the positions and type of atoms
xcoord = []; ycoord = []; zcoord = []; eltype = []; occfact = []; tempfact=[];


disp('Reading Atoms Positions...')
tic
for j = 1:1:size(rpdb.Model,2)
    if j >= 2
        shift = size(m,2);
    else
        shift = 0;
    end
    if (size(rpdb.Model,2) < 2)
        for ii=1:size(rpdb.Model.Atom,2)
            xcoord    = [xcoord rpdb.Model.Atom(1,ii).X];
            ycoord    = [ycoord rpdb.Model.Atom(1,ii).Y];
            zcoord    = [zcoord rpdb.Model.Atom(1,ii).Z];
            eltype{ii}= [rpdb.Model.Atom(1,ii).element];
            tempfact  = [tempfact; rpdb.Model.Atom(1,ii).tempFactor];
            occfact =[occfact; rpdb.Model.Atom(1,ii).occupancy];
            [~,~,m(ii)] = parametrizeScFac(eltype(ii));
        end
    else
        if isempty(rpdb.Model(j).Atom)
            for ii=1:size(rpdb.Model(j).HeterogenAtom,2)
                xcoord    = [xcoord rpdb.Model(j).HeterogenAtom(1,ii).X];
                ycoord    = [ycoord rpdb.Model(j).HeterogenAtom(1,ii).Y];
                zcoord    = [zcoord rpdb.Model(j).HeterogenAtom(1,ii).Z];
                eltype{shift+ii}= [rpdb.Model(j).HeterogenAtom(1,ii).element];
                tempfact  = [tempfact; rpdb.Model(j).HeterogenAtom(1,ii).tempFactor];
                occfact =[occfact; rpdb.Model(j).HeterogenAtom(1,ii).occupancy];
                [~,~,m(shift+ii)] = parametrizeScFac(eltype(shift+ii));
            end
        else
            for ii=1:size(rpdb.Model(j).Atom,2)
                xcoord    = [xcoord rpdb.Model(j).Atom(1,ii).X];
                ycoord    = [ycoord rpdb.Model(j).Atom(1,ii).Y];
                zcoord    = [zcoord rpdb.Model(j).Atom(1,ii).Z];
                eltype{shift+ii}= [rpdb.Model(j).Atom(1,ii).element];
                tempfact  = [tempfact; rpdb.Model(j).Atom(1,ii).tempFactor];
                occfact =[occfact; rpdb.Model(j).Atom(1,ii).occupancy];
                [~,~,m(shift+ii)] = parametrizeScFac(eltype(shift+ii));
            end
        end
    end
end
toc

% for ii=1:size(rpdb.Model.Atom,2)
%     xcoord    = [xcoord rpdb.Model.Atom(1,ii).X];
%     ycoord    = [ycoord rpdb.Model.Atom(1,ii).Y];
%     zcoord    = [zcoord rpdb.Model.Atom(1,ii).Z];
%     eltype{ii}= [rpdb.Model.Atom(1,ii).element];
%     tempfact  = [tempfact; rpdb.Model.Atom(1,ii).tempFactor];
%     occfact =[occfact; rpdb.Model.Atom(1,ii).occupancy];
%     [~,~,m(ii)] = parametrizeScFac(eltype(ii));
% end
    
    
xcoord = xcoord-min(xcoord);
ycoord = ycoord-min(ycoord);
Xc = xcoord*m'/sum(m);
Yc = ycoord*m'/sum(m);
xcoordCM = xcoord - Xc + abs(max(max(xcoord) - min(xcoord), max(ycoord) - min(ycoord)))/2;
ycoordCM = ycoord - Yc + abs(max(max(xcoord) - min(xcoord), max(ycoord) - min(ycoord)))/2;

xc0 = xcoordCM - min(xcoordCM) + extrapix;
yc0 = ycoordCM - min(ycoordCM) + extrapix;
szx0 = abs(max(xc0) - min(xc0)) + 2*extrapix;
szy0 = abs(max(yc0) - min(yc0)) + 2*extrapix;
mindist = max(szx0, szy0)*1/params2.acquis.pixsize*1e-10;

numpartExtra = params2.proc.partNum - params2.NumGenPart;
% limit the number of particles to the specified field of view 
if params2.proc.partNum> floor((params2.proc.N/(mindist)))^2
  fprintf('The maximum number of particles within this field of view is limited to %d\n', floor((params2.proc.N/(mindist)))^2)
  numpartExtra = floor((params2.proc.N/(mindist)))^2 - params2.NumGenPart;
end

outpot=[];
for ss = 1:numpartExtra  
    kkkk = tic;
    RAl    = [cos(alpha(ss)) sin(alpha(ss)); -sin(alpha(ss)) cos(alpha(ss))];
    RBet   = [cos(beta(ss)) sin(beta(ss)); -sin(beta(ss)) cos(beta(ss))];
    RGamma = [cos(gamma(ss)) sin(gamma(ss)); -sin(gamma(ss)) cos(gamma(ss))];
    x1y1 = RAl*[xcoordCM; ycoordCM];
    z2y2  = RBet*[zcoord; x1y1(2,:)];
    x3y3  = RGamma*[x1y1(1,:); z2y2(2,:)];
    xcoord3 = x3y3(1,:);
    ycoord3 = x3y3(2,:);
    zcoord3 = z2y2(1,:);

    xc = xcoord3 - min(xcoord3) + extrapix;
    yc = ycoord3 - min(ycoord3) + extrapix;
    zc = zcoord3 - min(zcoord3) + extrapix;
    % Define the volume of the protein
    szx = abs(max(xc) - min(xc)) + 2*extrapix;
    szy = abs(max(yc) - min(yc)) + 2*extrapix;
    szz = abs(max(zc) - min(zc)) + 2*extrapix;

    sz  = [szx szy szz]*oversam/voxel_size;
    ph   = phys_const;
    C   = 4*sqrt(pi)*ph.h^2/(ph.el*ph.me)*1e20; % =2132.8 A^2*V
    
    atompot = zeros(round(sz));
    
    fprintf('Calculating potential for %6d atoms of particle %3d\n ...',length(eltype),ss )
    for jj = 1:round(length(eltype))
%         if mod(jj, 5000)==0
%             fprintf('Calculating potential for the atom number %4d\n',jj)
%         end
        elem  = eltype(jj);
        occupancy = occfact(jj);
        [a,b,~] = parametrizeScFac(elem);
        b = b+tempfact(jj)+ 16*voxel_size^2;
        r2 = 0;
        b1 = zeros(1,5);
        for ll = 1:5
            b1(ll) = 4*pi^2/b(ll)*voxel_size^2;
            r2     = max(r2, 10/b1(ll)); % 10 corresponds to 4.5 sigma truncation
        end
        r   = sqrt(r2/3);
        xc1 = xc(jj)/voxel_size;
        yc1 = yc(jj)/voxel_size;
        zc1 = zc(jj)/voxel_size;
        rc  = [xc1, yc1, zc1];
        kmin = max(0, ceil((rc-r)));
        kmax = min(floor(sz)-1, floor((rc+r)));
        kmm  = max(kmax - kmin);

        x  = xc1 - [kmin(1): kmin(1)+kmm];
        y  = yc1 - [kmin(2): kmin(2)+kmm];
        z  = zc1 - [kmin(3): kmin(3)+kmm];
        x2 = x.^2;
        y2 = y.^2;
        z2 = z.^2;

        y2 = reshape(y2,[numel(y2) 1]);
        z2 = reshape(z2,[1 1 numel(z2)]);    
        newvalue = 0;
        for kk = 1:5
            tmp = a(kk)/b(kk)^(3/2) * C *repmat(exp(-b1(kk)*x2), [kmm+1 1 kmm+1]) .* repmat(exp(-b1(kk)*y2), [1 kmm+1 kmm+1]) .* repmat(exp(-b1(kk)*z2), [kmm+1 kmm+1 1]);        
            newvalue = newvalue + tmp;
        end
        newvalue=newvalue*occupancy;
        atompot(kmin(1):kmin(1)+kmm, kmin(2):kmin(2)+kmm, kmin(3):kmin(3)+kmm) = atompot(kmin(1):kmin(1)+kmm, kmin(2):kmin(2)+kmm, kmin(3):kmin(3)+kmm) + newvalue;
    end
    toc
    
    % fill in the empty space (with low potential values) with a contant value of vitreous ice which can be subtracted
    atompot1 = atompot - Vwat;
    atompot1(atompot1<0) = 0; 
    %it is neccesary to apply low-pass filter before eventual downsampling
    disp(' ')
    disp('Resampling...')
    tic
    if voxel_size < pixsz
        atompot1  = gaussf(mat2im(atompot1), sqrt((pixsz/voxel_size)^2-1), 'best');
    end
    atPotBlRspm = double(resample(atompot1, voxel_size/pixsz)); 
    toc
    
    disp(' ')
    disp('Adding Motion Blur...')
    tic
    % motion factor
    if params2.spec.motblur~=0
        atPotBlRspm = double(motionBlur(atPotBlRspm,params2));
    end
    toc
    
    %imaginary part
    if  params2.spec.imagpot == 3
    % fill in the empty space (up to the water potential value) with a contant value of vitreous ice which can be subtracted
    atompot2 = atompot - Vwat;
    atompot2(atompot2<0) = 0;  
    imagPot = params2.spec.potenampl*(atompot2==0) + params2.spec.proteinampl*(atompot2>0);
    imagPot = imagPot - params2.spec.potenampl;
    %it is neccesary to apply low-pass filter before eventual downsampling
    if voxel_size < pixsz
       imagPot  = gaussf(mat2im(imagPot), sqrt((pixsz/voxel_size)^2-1), 'best');
    end
    atPotBlRspmIm = double(resample(imagPot, voxel_size/pixsz));
    % motion factor
        if params2.spec.motblur~=0
            atPotBlRspmIm = double(motionBlur(atPotBlRspmIm,params2));
        end
        %atpotcmpx=atpot4blrspmWat+1i*atpot4blrspmWatIm;
    end  
  
    szPot = size(double(atPotBlRspm));
     if  params2.spec.imagpot == 3
       outputfilename = sprintf('%s_a%04d_Nx%i_Ny%i_Nz%i_Alp%3.1f_Bet%3.1f_Gam%3.1f_MF%3.1f_VoxSize%02.2fA_Volt%03dkV.raw',params2.spec.pdbin,ss+params2.NumGenPart,szPot, alphad(ss), betad(ss), gammad(ss),params2.spec.motblur, pixsz, params2.acquis.Voltage/1000);
    else
       outputfilename = sprintf('%s_a%04d_Nx%i_Ny%i_Nz%i_Alp%3.1f_Bet%3.1f_Gam%3.1f_MF%3.1f_VoxSize%02.2fA.raw',params2.spec.pdbin,ss+params2.NumGenPart,szPot, alphad(ss), betad(ss), gammad(ss),params2.spec.motblur, pixsz);
    end
   
    if ~exist([dir0 filesep 'Particles'])
        mkdir([dir0 filesep],'Particles')
    end
    
    OutFileName = [dir0 filesep 'Particles' filesep outputfilename];

    %disp(OutFileName);
    if wr
        fid = fopen(OutFileName,'w');
        if  params2.spec.imagpot ~= 3
            bla= fwrite(fid,double(atPotBlRspm), 'double');
        else
        % if the difference in mean free inelastic paths between protein and
        % vitreous is considered, the input 3D interaction potential is extended with
        % imaginary part
            bla= fwrite(fid,[double(atPotBlRspm); double(atPotBlRspmIm)], 'double');
        end
        %disp(bla)
        fclose(fid);
    end
    disp('Total time for calculating potential')
    toc(kkkk)
    disp('-----------------------------------------------------------------')
end