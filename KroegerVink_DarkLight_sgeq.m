%This code evaluates the defect chemistry of a mixed ionic-electronic
%conductor under dark and under light based on a 0 dimensional model
%(reaction limited, no transport included). The input parameters are
%relevant to halide perovskites such as methylammonium lead iodide (MAPI),
%with anti-Frenkel ionic (iodide) disorder, however it can be adapted to
%other material systems, by changing the relevant electronic and ionic
%properties. The surface exchange kinetics are not explicitly considered 
%(i.e. sg-eq condition is valid)

%How to use the code:
%- Select directory, folder name and calculation name
%- Change input parameters, such as Gext, Gamma_I_i, Gamma_I_v, Gamma_p_i,
%Gamma_n_v, recombination paramteres etc.
%- Make sure that the Equations_DarkLight_sgeq.m file is stored in the same
%folder as this script
%- Click on Run

%What the code gives as output:
%The code automatically plots relevant graphs with defect concentrations,
%QFLS and ionic chemical potential, etc. Additional plots can be defined.
%The code also saves txt files with the most important data and an info
%file with the input parameters used in the calculation

%--------------------------------------------------------------------------
% LICENSE
% Copyright (C) 2025  Davide Moia
% Max Planck Institute for Solid State Research
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%--------------------------------------------------------------------------
 
%Directory and Filename
Directory = 'Directory\';
newFolder = 'sg-eq';
[void] = mkdir(Directory,newFolder);
Filename = '\Filename';
NewDirectory = [Directory newFolder];

Include_dark = 1;       %Set to 0 if only solution under light is needed

rad_coeff = 1;       %Set to 0 to deactivate radiative recombination
nr_coeff = 1;        %Set to 0 to deactivate non-radiative recombination due to immobile defects (see parameters Nt, taun, taup, n1, p1 below)
aug_coeff = 1;       %Set to 0 to deactivate Auger recombination

i_coeff = 1;         %Set to 0 to if interstitials are redox inactive, 1 otherwise (within this version i_coeff should be set to 1, interstitials are by default contributing to the solid-gas excahnge)
v_coeff = 0;         %Set to 0 to if vacancies are redox inactive, 1 otherwise (within this version v_coeff can be 1 only if i_coeff is 1 too)

Rec_Model = [rad_coeff nr_coeff aug_coeff];  %Recombination model to use. The five terms can be 1 or 0 and they multiply the radiative, non-radiative (immobile trap), Auger, non-radiative (mobile trap), and iodine exchange terms
Redox_Model = [i_coeff v_coeff];             %Redox Model indicates if the reactions associated with the iodide intersitital (process i) and/or iodide vacancies (process v) should be used

color = {'r--' 'b--' 'g--' 'k--';'r.' 'b.' 'g.' 'k.'};

%Number of points N (values for pI2), ds is used later to downsample the data, so that all files have 401 lines
%By changing ds, the total number of calculations can be  varied. This can
%sometimes help with convergence.
ds = 4;                
N = ds*401;

%Values for normalization. These values can be varied to improve
%convergence
conc_ = 1e14;

%Temperature
T = 300;

%Constants
q = 1.6e-19;
Vth = 1.38e-23*T/q;

%Band gap, product of effective density of states for conduction and 
%valence band (assumed equal), intrinsic concentration of electrons and 
%holes
Eg = 1.63;
NcNv = 1e38;
ni = sqrt(NcNv)*exp(-Eg/2/Vth)/conc_;   %Electronic disorder mass action constant K_B = ni^2

%SRH-active trap parameters n1, p1, taun and taup (immobile defects) and
%density
n1 = ni;                                
p1 = ni^2/n1;
taun = 2e-6;
taup = 2e-6;
Nt = 1e18/conc_;                        

%Illumination paramters. The approximate generation rate is evaluated for
%MAPI, considering that a solar cell with 500 nm active layer yields about 
%23 mA/cm^2.  
SunsEq = 1e-3;                          %Illumination intensity in suns equivalent
Gext = SunsEq*0.023/q/5e-5/conc_; %External generation term, here considered constant in the material

krad = 1e-11*conc_;               %Radiative recombination constant 
Cn = 1e-28*conc_^2;               %Auger coefficients
Cp = 1e-28*conc_^2;

%Bimolecular generation-recombination of vacancies and interstitials is
%considered for the anti-Frenkel disorder
kion_rec = 1e-15*conc_;           %Bimolecular ionic defects recombination rate constant (kf_aF in the manuscript)
KaF = 1e32/conc_^2;                     %anti-Frenkel disorder mass action constant

PI2_IP = 1e-2;                %PI2 at the boundary between the I and P region (p = [Ii']). Estimate from conductivity measurements assuming mobility for n and p of 10 cm^2 V^-1 s^-1 and for ion 1e-9 cm^2 V^-1 s^-1
PI2_NI = PI2_IP/(KaF/ni^2)^2;           %PI2 at the boundary between the N and I region
PI2i = PI2_IP*ni^2/KaF;                 %PI2 corresponding to the intrinsic condition p = n

%Energy position above the valence band edge of the redox level associated 
%with iodide interstitials
E_Ii = 0.2979;
%Energy position below the conuction band edge of the redox level associated 
%with iodide vacancies (positive value required)
E_VI = 0.2979;

%Mass-action constants for the solid-gas exchange at the surface. The
%concentration of neutral defects (I_i^x at PI2_IP and V_I^x at PI2_NI) are
%calculated based on the energy position of the redox level associated with 
%each defect.
Iix_IP = KaF/sqrt(NcNv/conc_^2)*exp(E_Ii/Vth);
K_sg_i = sqrt(PI2_IP)/Iix_IP;            %K_sg_i = P(I2)^0.5/[I_i^x] defined at the boundary I P regions
VIx_NI = KaF/sqrt(NcNv/conc_^2)*exp(E_VI/Vth);
K_sg_v = sqrt(PI2_NI)*VIx_NI;           %K_sg_v = P(I2)^0.5*[V_I^x] defined at the boundary N I regions

%Process (i) is the redox with an iodide interstitial 
K_p_i = sqrt(PI2_IP)/KaF/K_sg_i;            %Pseudo mass action constant (at equilibrium), for the Ii'/Iix mediated by holes K_p_i = [Iix]/p/[Ii'] = pI2^0.5/Ki/p/[Ii']
K_n_i = K_p_i*ni^2;                     %Pseudo mass action constant (at equilibrium), for the Ii'/Iix mediated by electrons  K_n_i = [Iix]*n/[Ii'] = pI2^0.5*n/[Ii']/Ki = KI_p_i*ni^2 (from detailed balance)

%Process (v) is the redox with an iodide vacancy  
K_p_v = sqrt(PI2_IP)/K_sg_v;                %Pseudo mass action constant (at equilibrium), for the VIx/VI mediated by holes K_p_v = [VI.]/p/[VIx] = [VI.]pI2^0.5/Kv/p
K_n_v = K_p_v*ni^2;                     %Pseudo mass action constant (at equilibrium), for the VIx/VI mediated by electrons  K_n_v = [VI.]/p/[VIx] = [VI.]pI2^0.5/Kv/p = KI_h*ni^2 (from detailed balance)

vth = 1e7;                              %Electrons and holes thermal velocity

Gamma_p_i = 1e-2;                           %Parameter describing the ratio of the hole trapping rate due to redox with interstitials and the radiative rate evaluated at equilibrium and at P(I2)i Gamma_p_i = kf_p_i*KaF^0.5/(krad*ni)
CrossSection_Iim_p = Gamma_p_i*krad*ni/vth/KaF^0.5;         %Reference value for the cross section related with capture of a hole by a Ii'
CaptureCoeff_Iim_p = vth*CrossSection_Iim_p;

Gamma_n_v = 1e-2;                           %Parameter describing the ratio of the electron trapping rate due to redox with vacancies and the radiative rate evaluated at equilibrium and at P(I2)i Gamma_n_v = kf_p_v*Kv*P(I2)^-0.5/(krad*ni) 
CrossSection_VIp_n = Gamma_n_v*krad*ni/vth/KaF^0.5;         %Reference value for the cross section related with capture of an electron by a VI.
CaptureCoeff_VIp_n = vth*CrossSection_VIp_n;

%Gamma_I indicates the relative tendency of intestitials (i) or vacancies (v) to react with
%holes (Gamma_I>>1) or electrons (Gamma_I<<1)
Gamma_I_i = 1e0;                       %Ratio of the rates of Ii' oxidation due to hole capture vs electron release at P(I2)i at equilibrium Gamma_I_i = kf_p_i*ni/kf_n_i
Gamma_I_v = 1e0;                        %Ratio of the rates of VI. reduction due to hole release vs electron capture at P(I2)i at equilibrium Gamma_I_v = kf_p_v*ni/kf_n_v

kf_p_i = CaptureCoeff_Iim_p;            %Rate constant for Ii' oxidation via hole (process i)
kf_n_i = kf_p_i*ni*Gamma_I_i^-1;        %Rate constant for Ii' oxidation releasing an electron (process i)
kb_p_i = kf_p_i/K_p_i;                  %Rate constant for Iix reduction releasing a hole (process i)
kb_n_i = kf_n_i/K_n_i;                  %Rate constant for Iix reduction via an electron (process i)

kb_n_v = CaptureCoeff_VIp_n;            %Rate constant for VI. reduction via an electron (process v)
kb_p_v = kb_n_v*ni*Gamma_I_v;           %Rate constant for VI. reduction releasing a hole (process v)   
kf_n_v = kb_n_v*K_n_v;                  %Rate constant for VIx oxidation releasing an electron (process v). 
kf_p_v = kb_p_v*K_p_v;                  %Rate constant for VIx oxidation by a hole (process v)

%Partial pressure range. Decades_pI2_Nregion and Decades_pI2_Pregion indicate the number of
%decades in P(I2) for the P and the N regions, respectively
Decades_pI2_Nregion = 10;
Decades_pI2_Pregion = 10;
%Approximated boundary values for the intrinsic-N and intrinsic-P regions
pI2_intr_region_low = ni^4*(K_p_v*K_sg_v)^2/KaF^2;
pI2_intr_region_high = (K_p_v*K_sg_v)^2;
Log_pI2_min = log10(pI2_intr_region_low) - Decades_pI2_Nregion;
Log_pI2_max = log10(pI2_intr_region_high) + Decades_pI2_Pregion;
pI2 = logspace(Log_pI2_min, Log_pI2_max, N);
logstep = (Log_pI2_max - Log_pI2_min)/(N-1); 

%Calculate the solution under dark and then under light by solving the equations in the script
%Equations_DarkLight. The solution of the previous point is used as intial 
%guess for the next point

%Options for the solver 
%trust region dogleg, trust region, and Levenberg-Marquardt. 
%OPTIONS = optimoptions('fsolve','Algorithm','trust-region'),
options = optimoptions('fsolve','MaxIterations',1e18,'MaxFunctionEvaluations',1e18);

%Dark

%Initial guess for the solution at equilibrium at the lowest P(I2) calculated based on the Brouwer
%approximation. 
p0dark = ni^2/(KaF^0.5)*10^(-Decades_pI2_Nregion/4);
n0dark = KaF^0.5*10^(Decades_pI2_Nregion/4);
VI0dark = KaF^0.5*10^(Decades_pI2_Nregion/4);
Ii0dark = KaF^0.5*10^-(Decades_pI2_Nregion/4);

if Include_dark == 1
    
    xsol = [p0dark n0dark VI0dark Ii0dark];
    xsolall_dark = zeros(N,4);
    xguessall_dark = zeros(N,4);

    for i = 1:N 
%         i
        par = [pI2(i), kf_n_i, kf_p_i, kb_n_i, kb_p_i, kf_n_v, kf_p_v, kb_n_v, kb_p_v, KaF, kion_rec, ni, n1, p1, taun, taup, 0, krad, Cn, Cp, K_sg_i, K_sg_v];
        fun = @(x)Equations_DarkLight_sgeq(x, par, Rec_Model, Redox_Model);
        if i<3
            % Assuming a slope of 1/4 and -1/4 the initial guesses are
            % calculated for the first 2 data points
            x0 = [xsol(1)*10^((i-1)*logstep/4) xsol(2)*10^(-(i-1)*logstep/4) xsol(3)*10^(-(i-1)*logstep/4) xsol(4)*10^((i-1)*logstep/4)];
        else
            %For all other data points, the guess is linearly extrapolated from 
            %previous data points
            x0 = xsolall_dark(i-1,:).*((xsolall_dark(i-1,:)./xsolall_dark(i-2,:)).^(1))./conc_;
        end
        xsol = fsolve(fun, x0, options);
        xsolall_dark(i,:) = xsol.*conc_;
        xguessall_dark(i,:) = x0.*conc_;
    end

    %Plotting Dark solution and guesses
    Dash_type = 1;
    figure(1)
    loglog(pI2,xsolall_dark(:,1),color{Dash_type,1},pI2,xsolall_dark(:,2),color{Dash_type,2},pI2,xsolall_dark(:,3),color{Dash_type,3},pI2,xsolall_dark(:,4),color{Dash_type,4})
    %Optionally plotting the initial guesses for diagnostics in case of poor
    %convergence
%     hold on
%     loglog(pI2,xguessall_dark(:,1),'rx',pI2,xguessall_dark(:,2),'bx',pI2,xguessall_dark(:,3),'gx',pI2,xguessall_dark(:,4),'kx')
end

%Light

%Create the solution matrix. N rows (P(I2) points) and 4 columns for the variables [p n VI VMA]
xsolall = zeros(N,4);
xguessall = zeros(N,4);

%Paramater kn introduced to improve initial guess for n, if needed. kn is
%a guess of how much larger is the electron concentration under light than 
%in the dark for the lowest pI2
kn = 10;
p0 = max(p0dark,min([Gext/krad/n0dark/kn/Rec_Model(1),Gext*taup/Rec_Model(2),Gext/(n0dark*kn)^2/(Cn)/Rec_Model(3)]));
n0 = n0dark*kn;
VI0 = VI0dark;
Ii0 = Ii0dark;

% xsol = [max((ni^2/(KaF^0.5)*10^(-Decades_pI2_Nregion/4)),min([(Gext/krad/(KaF^0.5*10^(Decades_pI2_Nregion/4)))/Rec_Model(1),Gext*taup/Rec_Model(2),(Gext/(KaF^0.5*10^(Decades_pI2_Nregion/4))^2/(Cn))/Rec_Model(3)])) KaF^0.5*10^(Decades_pI2_Nregion/4) KaF^0.5*10^(Decades_pI2_Nregion/4) KaF^0.5*10^(-Decades_pI2_Nregion/4)];
xsol = [p0 n0 VI0 Ii0];

for i = 1:N 
%     i
    par = [pI2(i), kf_n_i, kf_p_i, kb_n_i, kb_p_i, kf_n_v, kf_p_v, kb_n_v, kb_p_v, KaF, kion_rec, ni, n1, p1, taun, taup, Gext, krad, Cn, Cp, K_sg_i, K_sg_v];
    fun = @(x)Equations_DarkLight_sgeq(x, par, Rec_Model, Redox_Model);
    %Initial guess based on extrapolation of previous two solution
    %datapoints
    if i<3
        % Assuming a slope of 1/4 and -1/4 the initial guesses are
        % calculated for the first 2 data points
        x0 = [xsol(1)*10^((i-1)*logstep/4) xsol(2)*10^(-(i-1)*logstep/4) xsol(3)*10^(-(i-1)*logstep/4) xsol(4)*10^((i-1)*logstep/4)];
    else
        %For all other data points, the guess is linearly extrapolated from 
        %previous data points
        x0 = xsolall(i-1,:).*((xsolall(i-1,:)./xsolall(i-2,:)).^(1))./conc_;
    end
    xsol = fsolve(fun, x0, options);
    xsolall(i,:) = xsol.*conc_;
    xguessall(i,:) = x0.*conc_;    
end

%Plotting the results

%Kroeger-Vink diagram
%Light
Dash_type = 2;
figure(1)
hold on
loglog(pI2,xsolall(:,1),color{Dash_type,1},pI2,xsolall(:,2),color{Dash_type,2},pI2,xsolall(:,3),color{Dash_type,3},pI2,xsolall(:,4),color{Dash_type,4});%pI2,pI2.^0.5/Ki,'ko',pI2,Kv./pI2.^0.5,'go')
%Optionally plotting the initial guesses for diagnostics in case of poor
%convergence
% loglog(pI2,xguessall(:,1),'r',pI2,xguessall(:,2),'b',pI2,xguessall(:,3),'g',pI2,xguessall(:,4),'k')
xlabel('P(I_2) (bar)')
ylabel('Defect concentration')

%QFLS as function of P(I2)
QFLS = Vth*log(xsolall(:,1).*(xsolall(:,2))/(ni*conc_)^2);
figure(2)
% hold on
semilogx(pI2,QFLS)
xlabel('P(I_2) (bar)')
ylabel('QFLS')

% Ionic defect pair chemical potential as function of P(I2)
Deltamuion_light = Vth*log(xsolall(:,3).*(xsolall(:,4))/(KaF*conc_^2));
figure(3)
% hold on
semilogx(pI2,Deltamuion_light,'c')
xlabel('P(I_2) (bar)')
ylabel('QFLSion')

% Quasi-chemical potential of iodine defined as mu_p - mu~_VI = Ap*log10(p/VI) 
% and -mu_n - mu~_VI = An*log10(p/VI) in the dilute limit
deltamu_p = Vth*log(xsolall(:,1)./xsolall_dark(:,1));
deltamu_n = Vth*log(xsolall(:,2)./xsolall_dark(:,2));
deltamu_VI = Vth*log(xsolall(:,3)./xsolall_dark(:,3));
deltamu_Ii = Vth*log(xsolall(:,4)./xsolall_dark(:,4));
% 
% figure(4)
% semilogx(pI2,deltamu_p,'r',pI2,deltamu_n,'b',pI2,deltamu_VI,'g',pI2,deltamu_Ii,'k')
% xlabel('P(I_2) (bar)')
% ylabel('deltamu')

% %Net recombination terms
Urad = (xsolall(:,2).*xsolall(:,1)-(conc_*ni)^2)*(krad/conc_);
USRH = (xsolall(:,2).*xsolall(:,1)-(conc_*ni)^2)./((taup)*(xsolall(:,2)+conc_*n1)+(taun)*(xsolall(:,1)+conc_*p1));
UAug = (xsolall(:,2).*xsolall(:,1)-(conc_*ni)^2).*(Cp/conc_^2*xsolall(:,1) + Cn/conc_^2*xsolall(:,2));

Upi = (kf_p_i/conc_*xsolall(:,1).*xsolall(:,4) - kb_p_i.*(pI2'.^0.5)/K_sg_i*conc_);
Upv = (kf_p_v*K_sg_v*xsolall(:,1).*(pI2'.^-0.5) - kb_p_v*xsolall(3));
Uni = kb_n_i*xsolall(:,2).*pI2'.^0.5/K_sg_i - kf_n_i*xsolall(:,4);
Unv = kb_n_v*xsolall(:,3).*xsolall(:,2)/conc_ - conc_*kf_n_v*K_sg_v*pI2'.^-0.5;
% figure(5)
% loglog(pI2,Urad,'k',pI2,USRH,'m',pI2,UAug,'g',pI2,Upi,'r--',pI2,Upv,'r.',pI2,Upv+Upi,'rs',pI2,Uni,'b--',pI2,Unv,'b.',pI2,Unv+Uni,'bx')
% xlabel('P(I_2) (bar)')
% ylabel('Net recombination rates')

GaF = Uni - Upi;
GaFx = -GaF;
% 
% figure(6)
% semilogx(pI2,GaF,'g',pI2,GaFx,'k')
% xlabel('P(I_2) (bar)')
% ylabel('Effective ionic defect generation rate')

% %Stoichiometry
% figure(7)
% loglog(pI2,xsolall_dark(:,4)-xsolall_dark(:,3),'k--',pI2,-xsolall_dark(:,4)+xsolall_dark(:,3),'k.-',pI2,xsolall(:,4)-xsolall(:,3),'r',pI2,-xsolall(:,4)+xsolall(:,3),'r.-')
% xlabel('P(I_2) (bar)')
% ylabel('Stoichiometry, delta')

%Saving the data. The file Info contains the relevant input parameters
fid = fopen([NewDirectory, Filename,'_sgeq_info.txt'],'wt');
Info = ['Suns = ' num2str(SunsEq) newline 'Gext = ' num2str(Gext*conc_) ' cm^-3 s^-1' newline 'Rec_Model = ' num2str(Rec_Model) newline 'Redox_Model = ' num2str(Redox_Model)...
    newline 'KaF = ' num2str(KaF*conc_^2) ' cm^-6'...
    newline 'Eg = ' num2str(Eg) ' eV' newline 'NcNv = ' num2str(NcNv*conc_) ' cm^-6' newline 'ni  = ' num2str(ni*conc_) ' cm^-3'...     
    newline 'kion_rec = ' num2str(kion_rec/conc_) ' cm^3 s^-1' newline 'PI2_IP = ' num2str(PI2_IP) ' bar' newline 'PI2_NI  = ' num2str(PI2_NI) ' bar'...     
    newline 'PI2i = ' num2str(PI2i) ' bar'...  
    newline 'K_p_i = ' num2str(K_p_i/conc_) ' cm^3' newline 'K_n_i = ' num2str(K_n_i*conc_) ' cm^-3'...
    newline 'K_p_v = ' num2str(K_p_v/conc_) ' cm^3' newline 'K_n_v = ' num2str(K_n_v*conc_) ' cm^-3'...
    newline 'vth = ' num2str(vth) ' cm/s'...
    newline 'Gamma_p_i = ' num2str(Gamma_p_i) newline 'Gamma_n_v = ' num2str(Gamma_n_v) ...
    newline 'Gamma_I_i = ' num2str(Gamma_I_i) newline 'Gamma_I_v = ' num2str(Gamma_I_v) ...
    newline 'kf_p_i = ' num2str(kf_p_i/conc_) ' cm^3 s^-1' newline 'kb_p_i = ' num2str(kf_p_i) ' s^-1'...
    newline 'kf_n_i = ' num2str(kf_n_i) ' s^-1' newline 'kb_n_i = ' num2str(kb_n_i/conc_) ' cm^3 s^-1'...
    newline 'kf_p_v = ' num2str(kf_p_v/conc_) ' cm^3 s^-1' newline 'kb_p_v = ' num2str(kf_p_v) ' s^-1'...
    newline 'kf_n_v = ' num2str(kf_n_v) ' s^-1' newline 'kb_n_v = ' num2str(kb_n_v/conc_) ' cm^3 s^-1'...
    newline 'K_sg_i = ' num2str(K_sg_i/conc_) ' bar^1/2 cm^3' newline 'K_sg_v  = ' num2str(K_sg_v*conc_) ' bar^1/2 cm^-3' ...
    newline 'taun = ' num2str(taun) ' s' newline 'taup = ' num2str(taup) ' s'...
    newline 'n1 = ' num2str(n1*conc_) ' cm^-3' newline 'p1  = ' num2str(p1*conc_) ' cm^-3'... 
    newline 'krad = ' num2str(krad/conc_) ' cm^3 s^-1' newline 'Cn = ' num2str(Cn/conc_^2) ' cm^6 s^-1' newline 'Cp = ' num2str(Cp/conc_^2) ' cm^6 s^-1'...
    newline 'Decades_pI2_Nregion = ' num2str(Decades_pI2_Nregion) newline 'Decades_pI2_Pregion = ' num2str(Decades_pI2_Pregion) ...
    newline 'ds = ' num2str(ds) newline 'N = ' num2str(N)];
fprintf(fid, Info);
fclose(fid);

filename = [NewDirectory, Filename,'_sgeq_DefConc_light.txt'];
fid = fopen(filename, 'w');
fprintf(fid, '%s,%s,%s,%s,%s\n','log10(pI2/bar)','log10(p/cm^-3)','log10(n/cm^-3)','log10([V_I^.]/cm^-3)','log10([I_i^m]/cm^-3)');
fclose(fid);
dlmwrite(filename,[downsample(log10(pI2'),ds) downsample(log10(xsolall),ds)],'-append');

filename = [NewDirectory, Filename,'_sgeq_Deltamu.txt'];
fid = fopen(filename, 'w');
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s\n','log10(pI2/bar)','deltamu_p/eV','deltamu_n/eV','deltamu_V_I^./eV','deltamu_I_i^m/eV','deltamu_Ip/eV','deltamu_In/eV');
fclose(fid);
dlmwrite(filename,[downsample(log10(pI2'),ds) downsample(deltamu_p,ds) downsample(deltamu_n,ds) downsample(deltamu_VI,ds) downsample(deltamu_Ii,ds) downsample(deltamu_p-deltamu_VI,ds) downsample(-deltamu_n-deltamu_VI,ds)],'-append');

filename = [NewDirectory, Filename,'_sgeq_QFLS.txt'];
fid = fopen(filename, 'w');
fprintf(fid, '%s,%s\n','log10(pI2/bar)','QFLS/eV');
fclose(fid);
dlmwrite(filename,[downsample(log10(pI2'),ds) downsample(QFLS,ds)],'-append');

filename = [NewDirectory, Filename,'_sgeq_Deltamu_ion.txt'];
fid = fopen(filename, 'w');
fprintf(fid, '%s,%s\n','log10(pI2/bar)','QFLS/eV');
fclose(fid);
dlmwrite(filename,[downsample(log10(pI2'),ds) downsample(Deltamuion_light,ds)],'-append');

filename = [NewDirectory, Filename,'_sgeq_Recombination.txt'];
fid = fopen(filename, 'w');
fprintf(fid, '%s,%s,%s,%s,%s,%s\n','log10(pI2/bar)','Urad/cm^-3 s^-1','USRH/cm^-3 s^-1','UAug/cm^-3 s^-1','Upi/cm^-3 s^-1','Upv/cm^-3 s^-1');
fclose(fid);
dlmwrite(filename,[downsample(log10(pI2'),ds) downsample(log10(Urad),ds) downsample(log10(USRH),ds) downsample(log10(UAug),ds) i_coeff*downsample(log10(Upi),ds) v_coeff*downsample(log10(Upv),ds)],'-append');

filename = [NewDirectory, Filename,'_sgeq_EffectiveIonicGeneration.txt'];
fid = fopen(filename, 'w');
fprintf(fid, '%s,%s,%s\n','log10(pI2/bar)','G_aF/cm^-3 s^-1','G_aF^x/cm^-3 s^-1');
fclose(fid);
dlmwrite(filename,[downsample(log10(pI2'),ds) downsample(GaF,ds) downsample(GaFx,ds)],'-append');

if Include_dark == 1
    filename = [NewDirectory, Filename,'_sgeq_DefConc_darkequil.txt'];
    fid = fopen(filename, 'w');
    fprintf(fid, '%s,%s,%s,%s,%s\n','log10(pI2/bar)','log10(p/cm^-3)','log10(n/cm^-3)','log10([V_I^.]/cm^-3)','log10([I_i^m]/cm^-3)');
    fclose(fid);
    dlmwrite(filename,[downsample(log10(pI2'),ds) downsample(log10(xsolall_dark),ds)],'-append');
end
