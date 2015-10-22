function [noiseless_tilt_series, ctf_out]= DsimTEM2_ctf(psi_exit,params2)
%  CTF with df, ast, envelopes and optionally phase plate


% prealocate memory for the stack (series)
nTiltAngles = length(params2.acquis.tilt);

if strcmp(params2.seriesout,'tilt')
    Nseries= nTiltAngles;
elseif strcmp(params2.seriesout,'defocus')
    Nseries= length(params2.acquis.df);
elseif strcmp(params2.seriesout,'dose')
    Nseries= length(params2.acquis.dose_on_sample);
else
    Nseries=1;
end
% extprojstack  = newim(params2.proc.N,params2.proc.N,Nseries, 'complex');
btot_i        = newim(params2.proc.N,params2.proc.N,Nseries);
% series        = newim(params2.proc.N,params2.proc.N,Nseries);
psi_exit=double(psi_exit);   
btot_i=double(btot_i);

for jjj= 1:Nseries
    if strcmp(params2.seriesout,'defocus')
        params2.acquis.df_run=params2.acquis.df(jjj);
    else
        params2.acquis.df_run=params2.acquis.df(1);
    end
    if ~mod(jjj,5)||~mod(jjj,Nseries)
        fprintf(['Calculate the CTF for the ' params2.seriesout sprintf(' series. Image number %3d of %3d\n',  jjj, Nseries)]);
    end
  
  [ctf] = simulateCTF(params2);
  if params2.disp.ctf 
    ctf_out(:,:,jjj)=double(ctf);
  end
  btot = ctf*dip_fouriertransform(dip_image(psi_exit(:,:,jjj)),'forward',[1 1 ]); %0]);
  btot_i(:,:,jjj) = double(abs(dip_fouriertransform(btot,'inverse',[1 1])).^2); % intensity in the image without camera influence              
end
noiseless_tilt_series = dip_image(btot_i); 