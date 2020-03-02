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
% Nz = size(volume,2);   %added in 29/10/2019
nangles = numel(angles);

% allocate output memory
stack = newim(Nx,nangles,Ny);
stack = complex(stack,stack);

% qy = yy(1,Ny,'freq'); %orignial code

% FT in y direction
volume = dip_fouriertransform(volume,'forward',[0 0 1]);
volume = double(volume);
stack = double(stack);
% qy    = double(qy); %original code


% try to iterate over z
% if angles == 0                   %added on 29/19/2019
%     
%     qz = yy(1,Nz,'freq');
%     qz    = double(qz);
%     
%     parfor k = 1:Nz    
%     if ~mod(Nz-k,50)|| ~mod(k,Nz)
%         fprintf('Calculate Radon for the slice number %3d of %3d\n',  Nz-k, Nz)
%     end
%     slice = squeeze(volume(k,:,:));
%     
%     sinogram = FreqDomRadon_Ewald_2d_slice(dip_image(slice),qz(k),angles,params);
%     
%     stack(k,:,:) = double(sinogram);
%     end
% end  
%     
% else

% iterate over y    
    qy = yy(1,Ny,'freq');
    qy = double(qy);
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
