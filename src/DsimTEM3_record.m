function series = DsimTEM3_record(noiseless_tilt_series,params2)
%  --------------------------------- Camera influence  

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

series        = newim(params2.proc.N,params2.proc.N,Nseries);       
noiseless_tilt_series=double(noiseless_tilt_series);
series = double(series);

for iii= 1:Nseries
    if strcmp(params2.seriesout,'dose')
        params2.influx=params2.acquis.dose(iii);
    else
        params2.influx=params2.acquis.dose(1);
    end
    if ~mod(iii,5)||~mod(iii,Nseries)
        fprintf(['Calculate the DQE for the ' params2.seriesout sprintf(' series. Image number %3d of %3d\n', iii, Nseries)]);
    end
IntenDetect = DetectSim(squeeze(dip_image(noiseless_tilt_series(:,:,iii))), params2); 
series(:,:,iii) = double(IntenDetect);
end
series = dip_image(series);