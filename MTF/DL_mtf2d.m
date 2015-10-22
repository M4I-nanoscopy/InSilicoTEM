%%

dirMTF = [pwd, filesep];
f2mtf1d = myload([dirMTF,'MTF.mat']);
f2dqe1d = myload([dirMTF,'DQE.mat']);
figure(1)
subplot(121), plot(f2mtf1d(:,1), f2mtf1d(:,2), '-r');
subplot(122), plot(f2dqe1d(:,1), f2dqe1d(:,2), '-b');
%%

mtf1dx = f2mtf1d(:,1);
mtf1dy = f2mtf1d(:,2);
dqe1dx = f2dqe1d(:,1);
dqe1dy = f2dqe1d(:,2);

%% 
mtf1d = f2mtf1d(:,2);
mtf1d(mtf1d>1) = 1;
% resmple;
mtf1d128= resample(mtf1d, 128/length(mtf1d) );
mtf2d = from1d2d(mtf1d128,256);
mtf2d(mtf2d<=0) = mtf1d(end);
mtf = mtf2d;
% mtf1d = resample(EMGcrop,128/size(EMGcrop,2));
%
t128 = double(mtf1d128);
t2d = double(mtf2d);
%%
dqe1d = f2dqe1d(:,2);
dqe1d128 = resample(dqe1d,128/length(dqe1d));
dqe2d = from1d2d(dqe1d128,256);
dqe = dqe2d + dqe1d(end)*(dqe2d<=0);
tq2d = double(dqe);
%%
figure(2);
subplot(121), plot(mtf1d128,'-r');
subplot(122), plot(dqe1d128,'-b');

%%
dirSave = dirMTF;
save([dirSave,'MTF_FalconII_200.mat'],'mtf');
save([dirSave,'DQE_FalconII_200.mat'],'dqe');
%%
% mu = 0;
% L = 0.8;  % smaller the L closer to the exponential (here we will keep it constant)
% %s= 1;   % larger the sigma broader the distribution
% s = 120e3/params.acquis.Voltage; % lower the voltage, larger the sigma, better the MTF, and closer to Gaussian distibution
% x = -6:0.1:6;
% EMG = L/2*exp(L/2*(2*mu+L*s^2-2.*x)).*erfc((mu+L*s^2-x)/(sqrt(2)*s));
% [maxval,maxpos] = max(EMG);
% EMGcrop = EMG(maxpos:end)/maxval;
% mtf1d = resample(EMGcrop,128/size(EMGcrop,2));
% mtf = from1d2d(mtf1d,256);
% SimMTF = mtf+EMGcrop(end)*(mtf<=0);