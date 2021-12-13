function [circles, Nrparticles] = GetParticleCoordinate(mindist, params)

randdimension = params.spec.partdistribution;
pix = params.proc.N; % in A
thickness = (params.spec.thick)/params.acquis.pixsize; % in A
ParticleDensity = params.spec.partdensity;
tiltmax = max(params.acquis.tilt);

tic
disp('Generating random particle positions...')
if strcmp(randdimension, '2D')
    circles = Randomposition(pix,mindist,tiltmax);
    [row column] = size(circles);
    Nrparticles = row;
    
elseif strcmp(randdimension, '3D')
     if strcmp(params.inter.type,'wpoa')
%     circles = RandomPosition3D([pix, pix, thickness],mindist/2, ParticleDensity);
    circles = RandomPosition3D([floor((pix+mindist)/cos(tiltmax)),pix+mindist, thickness],mindist/2, ParticleDensity);
    [row column] = size(circles);
    Nrparticles = row;
     else
    circles = RandomPosition3D([pix+mindist, floor((pix+mindist)/cos(tiltmax)), thickness],mindist/2, ParticleDensity);
    [row column] = size(circles);
    Nrparticles = row;
     end
end
toc
end

    