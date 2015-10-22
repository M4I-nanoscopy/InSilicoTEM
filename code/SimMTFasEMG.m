function SimMTF = SimMTFasEMG(params)
%SimMTFasEMG Simulates the modulation transfer function (MTF) as
% exponentially modified Gaussian (EMG). For lower voltages, the MTF looks
% more like Gaussian, for higher voltages as exponential function
% SYNOPSIS:
% SimMTF = SimMTFasEMG(params)
%
% PARAMETERS:
%   params: structure containing various input paramters (e.g. Voltage)
%
% OUTPUT:
%  SimMTF: simulated MTF

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

mu = 0;
L = 0.8;  % smaller the L closer to the exponential (here we will keep it constant)
%s= 1;   % larger the sigma broader the distribution
s = 120e3/params.acquis.Voltage; % lower the voltage, larger the sigma, better the MTF, and closer to Gaussian distibution
x = -6:0.1:6;
EMG = L/2*exp(L/2*(2*mu+L*s^2-2.*x)).*erfc((mu+L*s^2-x)/(sqrt(2)*s));
[maxval,maxpos] = max(EMG);
EMGcrop = EMG(maxpos:end)/maxval;
mtf1d = resample(EMGcrop,128/size(EMGcrop,2));
mtf = from1d2d(mtf1d,256);
SimMTF = mtf+EMGcrop(end)*(mtf<=0);
%     plot(SimMTF)