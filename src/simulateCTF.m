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
% 
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
if df_run == 0 && params2.acquis.ast == 0
    inratioqqs=1;
    inratioqlq=1; 
end
    xdot = x*cos(params2.acquis.astangle) - y*sin(params2.acquis.astangle);
    ydot = x*sin(params2.acquis.astangle) + y*cos(params2.acquis.astangle);
    q = q + sqrt((xdot/inratioqlq ).^2+(ydot*inratioqqs).^2)*q_true_pix;
    qsym = qsym + sqrt((xdot).^2+(ydot).^2)*q_true_pix;
    c1= 2*pi*( 0.25*params2.mic.Cs*params2.inter.lambda^3*qsym.^4 - 0.5*df_run*params2.inter.lambda*q.^2);
    c = cos(c1)-1i*sin(c1);
    if isfield(c, 'data')
       bla = c.data;
       c = dip_image(bla);
    end
    % damping envelopes:
    
    % Kc is the envelope function due to chromatic aberration of the electron gun
    % Partial temporal coherence: Ec(u) = exp[-0.5*(pi*lamda*sigma)^2*u^4]
    % where sigma = Cc*sqrt[4*(deltaI/I)^2+(deltaE/voltage)^2+(deltaV/voltage)^2]
    H = params2.mic.C_c* params2.mic.deltaE/params2.acquis.Voltage;
    nomc = pi*params2.inter.lambda.*qsym.^2*H;
    denomc = 4*sqrt(log(2));
%    denomc = sqrt(2);
    Kc = exp(-(nomc/denomc).^2);
    if isfield(Kc, 'data')
       bla=Kc.data;
       Kc=dip_image(bla);
    end
    
    % Ks is the envelope function due to different direction of the electron coming out from the gun
    % Partial spatial coherence: Es(u) = exp[-(pi*alpha/lamda)^2*(Cs*lamda^3*u^3+lamda*u)^2]
    nums = (pi*params2.mic.Cs*params2.inter.lambda^2.*qsym.^3-pi*df_run.*q).^2*params2.mic.a_i^2;   
    Ks = exp(-nums/log(2));
%    Ks = exp(-nums);
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
    ctf = K.*c.*A; %total CTF
    
    % optional Spiral phase plate
    if params2.mic.SPPflag
%         central_blocked_PP = ones(params2.proc.N,params2.proc.N);
%         central_blocked_PP([(params2.proc.N/2-params2.proc.N/8):(params2.proc.N/2+1+params2.proc.N/8)] ...
%             , [(params2.proc.N/2-params2.proc.N/8):(params2.proc.N/2+1+params2.proc.N/8)])= 0;
%         Image_com=dip_image(central_blocked_PP,'complex');

        Image_com=dip_image(ones(params2.proc.N,params2.proc.N),'complex');
        
        % For spiral phase plate
        x_ic = xx(Image_com);
        y_ic = yy(Image_com);
        phi_ic = atan2(y_ic,x_ic);   
        SPhPlate = exp(1i*params2.mic.SPP_Phase*phi_ic).*Image_com; 
        
        SPhPlate(q < params2.mic.SPP_qcuton)= 1; %block the central beam = 0.1
        SPhPlate(qsym > qmax)= 0;
        SPhPlateIm = gaussf(imag(SPhPlate),1);
        SPhPlateRe = gaussf(real(SPhPlate),1);
        SPhPlate = SPhPlateRe+1i*SPhPlateIm;
        ctf = SPhPlate*ctf;
%       ctf = ctf*CTF_Envelope;
    end
     
% Optional Zernike Phase Plate 

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

% plot ctf (Spiral shape)

if params2.disp.ctf
    if params2.mic.SPPflag
        dipshow(ctf);
        ctf_array = dip_array(ctf);
   
        [r_ctf,c_ctf] = size(ctf_array);
        center_row = floor(r_ctf / 2);
        center_col = floor(c_ctf / 2);

        [x_ctf,y_ctf] = meshgrid((1:c_ctf) - center_col, (1:r_ctf) - center_row);

        [phi_ctf, rad_ctf] = cart2pol(x_ctf, y_ctf);

        r_ctf_range = linspace(0, params2.proc.N/2, params2.proc.N/2);
        phi_ctf_range = linspace(-pi, pi, params2.proc.N);

        [T, R] = meshgrid(phi_ctf_range, r_ctf_range);

        ctf_polar = griddata(phi_ctf, rad_ctf, real(ctf_array), T, R);
    
        plot([1:params2.proc.N/2]*q_true_pix/1e9, ctf_polar(:,params2.proc.N/2));
        hold on
        plot([1:params2.proc.N/2]*q_true_pix/1e9, ctf_polar(:,params2.proc.N/4*3), '-*');
        hold on
        plot([1:params2.proc.N/2]*q_true_pix/1e9, ctf_polar(:,params2.proc.N/4), '-o');

        ylabel('Contrast transfer function');
        xlabel(' Frequency [nm^-1]');
        ylim([-1,1]);
        title('1d CTF')
        legend('0.5*pi','0', 'pi')
        hold off
    end
end





%     Original code ctf plot: radialmean           

ctf1d = radialmean(imag(ctf));
    if params2.disp.ctf || params2.mic.SPPflag == 1
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
end
