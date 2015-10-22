function stack = FreqDomRadon_Ewald_2d(volume,angles,params)
%FREQDOMRADON_EWALD_2D   
%
%
% SYNOPSIS:
%  stack = FREQDOMRADON_EWALD_2D(volume,angles,params);
% 
% PARAMETERS:
%  volume:     Input volume (x,z,y)
%  angles:     Angles at which to compute the projections
%  params:     set of microscope and acquisition parameters
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
    error('stack should be 3D');
end

Nx = size(volume,1);
Ny = size(volume,3);
nangles = numel(angles);

% allocate output memory
stack = newim(Nx,nangles,Ny);
stack = complex(stack,stack);

qy = yy(1,Ny,'freq');

% FT in y direction
volume = dip_fouriertransform(volume,'forward',[0 0 1]);
volume = double(volume);
stack = double(stack);
qy    = double(qy);

% iterate over y
parfor k = 1:Ny    
    if ~mod(Ny-k,50)|| ~mod(k,Ny)
        fprintf('Calculate Radon for the slice number %3d of %3d\n',  Ny-k, Ny)
    end
    slice = squeeze(volume(:,:,k));
    
    sinogram = FreqDomRadon_Ewald_2d_slice(dip_image(slice),qy(k),angles,params);
    
    stack(:,:,k) = double(sinogram);
end

% IFT in y direction
stack = dip_fouriertransform(stack,'inverse',[0 0 1]);
