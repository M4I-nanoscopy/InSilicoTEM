function [x y z alpha beta gamma] = PartList(params2)
%PartList Specifies the translation and/or orientation of the particles
% within the volume if params.proc.geom = 1. In case it is 0 random
% orientation will be assumed. 
% SYNOPSIS:
% [x y z alpha beta gamma] = PartList(params2)
%
% PARAMETERS:
% params2: structure containing various input paramters (e.g. # of particles)
%
% OUTPUT:
%      x: translation in x direction [in pixles] from the image center
%      y: translation in y direction [in pixles] from the image center
%      z: translation in z direction [in pixles] from the image center
%  alpha: the first Euler angle
%   beta: the second Euler angle
%  gamma: the third Euler angle

n = params2.proc.partNum; 

% matrix with translations and rotations of particles
% the number of rows correspond to the number of particles
% the columns represent translation in x y z and three Euler rotations 
p  =  [   0   0    0      0   0   0;
          0   0    0      0   0   0;
          0   0    0      0   0   0;
          0   0    0      0   0   0;
          0   0    0      0   0   0;
          0   0    0      0   0   0;
          0   0    0      0   0   0;
          0   0    0      0   0   0;
          0   0    0      0   0   0;
          0   0    0      0   0   0];
     
 if size(p,1)<n
    warning(['At the moment the orientations/positions of only ' sprintf('%4d', size(p,1)) ' particles are specified. Please expand the PartList.m'])
 end
x = []; y = []; z = []; alpha = []; beta  = []; gamma = [];
for ii = 1:n
    x = [ x p(ii,1)];
    y = [ y p(ii,2)];
    z = [ z p(ii,3)];
    alpha = [alpha p(ii,4)];
    beta  = [beta p(ii,5)];
    gamma = [gamma p(ii,6)];
end



