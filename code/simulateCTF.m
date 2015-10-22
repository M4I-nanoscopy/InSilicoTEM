function [ctf] = simulateCTF(params2)
%simulateCTF simulates the effective contrast transfer function (CTF) by
% taking into account microscope
% aberrations and enevelopes due to uncoherency of the source
% SYNOPSIS:
% [ctf] = simulateCTF(params2)
%
% PARAMETERS:
% params2: structure containing various input physical and processing parameters
%
% OUTPUT:
%    ctf: effective contrast transfer function

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic


q_true_pix = 1/(params2.acquis.pixsize*params2.proc.N);
q          = newim([params2.proc.N params2.proc.N],'dfloat');
qsym       = newim([params2.proc.N params2.proc.N],'dfloat');
df_run     = params2.acquis.df_run;
x          = xx(q);
y          = yy(q);
dfsh       = df_run+params2.acquis.ast;
dfln       = df_run-params2.acquis.ast;
inratioqqs = sqrt(dfsh/df_run);
inratioqlq = sqrt(df_run/dfln);
if df_run==0 && params2.acquis.ast==0
    inratioqqs=1;
    inratioqlq=1;
end
xdot = x*cos(params2.acquis.astangle) - y*sin(params2.acquis.astangle);
ydot = x*sin(params2.acquis.astangle) + y*cos(params2.acquis.astangle);
q = q + sqrt((xdot/inratioqlq ).^2+(ydot*inratioqqs).^2)*q_true_pix;
qsym = qsym + sqrt((xdot).^2+(ydot).^2)*q_true_pix;
c1= 2*pi*( 0.25*params2.mic.Cs*params2.inter.lambda^3*qsym.^4 - 0.5*df_run*params2.inter.lambda*q.^2);
c = cos(c1)-1i*sin(c1); % ctf function
if isfield(c, 'data')
    bla = c.data;
    c = dip_image(bla);
end
% damping envelopes:
H = params2.mic.C_c* params2.mic.deltaE/params2.acquis.Voltage;
nomc = pi*params2.inter.lambda.*qsym.^2*H;
denomc = 4*sqrt(log(2));
Kc = exp(-(nomc/denomc).^2);
if isfield(Kc, 'data')
    bla=Kc.data;
    Kc=dip_image(bla);
end
nums = (pi*params2.mic.Cs*params2.inter.lambda^2.*qsym.^3-pi*df_run.*q).^2*params2.mic.a_i^2;
Ks = exp(-nums/log(2));
if isfield(Ks, 'data')
    bla=Ks.data;
    Ks=dip_image(bla);
end
K = Kc.*Ks;
% aperture function
A = dip_image(ones(params2.proc.N,params2.proc.N),'complex');
qmax= 2*pi*params2.mic.diam_obj/(params2.inter.lambda*params2.mic.foc);
%qmax=10*1e10;
A(q >qmax)=0;
A = gaussf(double(A),3); % in another code this is modeled as erf function and it is related to the focal distance (add that to this code!)
ctf = K.*c.*A; % total CTF
% optional phase plate
if params2.mic.PPflag
    PhPlate = exp(1i*params2.mic.PP_Phase)*dip_image(ones(params2.proc.N,params2.proc.N),'complex');
    PhPlate(q < params2.mic.PP_qcuton)=1;
    %PhPlate(rr (PhPlate) < 5)=1;
    PhPlate(qsym > qmax)=0;
    PhPlateIm = gaussf(imag(PhPlate),1);
    PhPlateRe = gaussf(real(PhPlate),1);
    PhPlate = PhPlateRe+1i*PhPlateIm;
    ctf = PhPlate*ctf;
end

% plot ctf
ctf1d = radialmean(imag(ctf));
if params2.disp.ctf
    plot([1:length(ctf1d)]*q_true_pix/1e9, ctf1d);
    hold on
    plot([1:length(ctf1d)]*q_true_pix/1e9, radialmean(K), 'g', [1:length(ctf1d)]*q_true_pix/1e9, radialmean(-K), 'g');
    ylabel('Contrast transfer function');
    xlabel(' Frequency [nm^-1]');
    %xlim( [0 4]);
    ylim( [-1 1]);
    title ('1D CTF')
    hold on
    plot([1:length(ctf1d)]*q_true_pix/1e9, 0, 'r');
    hold off
end

