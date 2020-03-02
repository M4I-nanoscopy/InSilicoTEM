function sinogram = FreqDomRadon_1d(slice,angles)
%FREQDOMRADON_1D   Compute Radon transform in Frequency domain
%
%
% SYNOPSIS:
%  sinogram = FREQDOMRADON_1D(slice,angles);
% 
% PARAMETERS:
%  slice:   Input slice (x,z)
%  angles:  Angles at which to compute the projections

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

% make sure input is dip_image
slice = dip_image(slice);

N = size(slice,1);
Nz = size(slice,2);

% create set of frequency points
qx = double(xx(N,1,'freq'));
qz = zeros(size(qx));

q = zeros(2,N*length(angles));

for t = 1:length(angles);
    R = [cos(angles(t)) sin(angles(t));-sin(angles(t)) cos(angles(t))];
    
    q(:,(t-1)*N+1:t*N) = R*[qx; qz];
end
q = dip_image(q);

% compute NUFFT
sinogram = nufft_type2(slice,im2array(q));
sinogram = reshape(sinogram,[N length(angles)]);

%inverse fourier transform
sinogram = real(dip_fouriertransform(sinogram,'inverse',[1 0]));
sinogram = 1/(2*N*Nz^(1/2))*sinogram;