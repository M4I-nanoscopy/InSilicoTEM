function stack = FreqDomRadon_2d(volume,angles)
%FREQDOMRADON_2D   Compute Radon transform in Frequency domain
%
%
% SYNOPSIS:
%  stack = FREQDOMRADON_2D(volume,angles);
% 
% PARAMETERS:
%  volume:   Input volume (x,z,y)
%  angles:  Angles at which to compute the projections
% 
% OUTPUT:
%  stack: 

% The Quantitative Imaging Group of the TU Delft has developed and is the 
% owner of the image processing routines distributed in this package, 
% hereafter called SOFTWARE.
% 
% The SOFTWARE is free for non-commercial use by students and staff in 
% universities or non-profit research institutes.
% 
% Redistribution of SOFTWARE or parts thereof in any form is not permitted.
% 
% This SOFTWARE is distributed in the hope that it will be useful, but 
% without any warranty; without even the implied warranty of 
% merchantability or fitness for a particular purpose.

volume = dip_image(volume);

if ndims(volume) ~= 3
    error('volume should be 3D');
end

% Nx = max(size(volume,1),size(volume,2));
Nx = size(volume,1);
Ny = size(volume,3);
nangles = numel(angles);

% allocate output memory
stack = newim(Nx,nangles,Ny);

% iterate over y
volume = double(volume);
stack = double(stack);
parfor k = 1:Ny
    if ~mod(Ny-k,50)|| ~mod(k,Ny)
        fprintf('Calculate Radon for the slice number %3d of %3d\n',  Ny-k, Ny)
    end
    slice = squeeze(volume(:,:,k));
    
    sinogram = FreqDomRadon_1d(dip_image(slice),angles);
    
    stack(:,:,k) = double(sinogram);
end
stack = dip_image(stack);


