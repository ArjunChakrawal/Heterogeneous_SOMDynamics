clc;close all;clearvars;

ksmm=0.018;
kb=0.00028;  %h-1
Y=0.31;
km=25;
rho= 1.65; %g/cm3 soil bulk density 
rho_wood=1.1;%g/cm3  max carbon density in form of wood 
nx=100;ny=100; %One pixel = 25micron

Vol_mic = 50*50*50; %um^3
Vmic_cc=Vol_mic*1e-12; %cm3
Vdomain_cc=Vol_mic*nx*ny*1e-12 ; % cm^3

maxC_ms=(rho_wood*1e-12)*Vol_mic*1e15;   %fg C
maxCD=maxC_ms*nx*ny; % fg C / vol domain 
gC_CC_maxCD=maxCD*1e-15/Vdomain_cc;% gC/ cc
gC_gS_maxCD=gC_CC_maxCD/rho;

max_allowed_ms=maxC_ms.*0.5;
total_allowedC= 1.5e7*(nx*ny);
avgC=total_allowedC/(nx*ny);

soil_domain=(rho*Vdomain_cc);
C_domain=total_allowedC*1e-15;
fac=C_domain/soil_domain;

%% ksmult estimation
load('Ph_mult_Csn31.mat');
load('Ph_mult_Cbn31.mat');
cs0=mean2(spCs)*1e-12./(rho*Vmic_cc);cb0=mean2(spCb)*1e-12./(rho*Vmic_cc);co20=0;
Imic=cs0/200000;

options = odeset('Stats','off','AbsTol',1e-6,'RelTol',1e-6);
T=24*1000; nstep=500; t=linspace(0,T,nstep);
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(t)).*Imic;
f=@(tt,c)[-ksmm *c(1)*c(2)/(km+c(1)) + kb*c(2) + interp1(t, St,tt);...
    Y*ksmm *c(1)*c(2)/(km+c(1)) - kb*c(2); (1-Y)*ksmm*c(1)*c(2)/(km+c(1))];
[~, Ch]=ode45(f,t,c(1,:), options); %homogeneous solution 

ks=(ksmm.*Y - kb)./(Y.*km);

p0=ks;
options = optimset( 'MaxFunEvals', 1000,'MaxIter', 1000, 'TolFun', 1e-10, 'TolX',1e-10);
[p_estimate,fval,exitflag,output]  = fminsearch(@(p)odefit(t,Ch,p,ksmm,km,kb,'mult',St),p0,options);

%%
my=(ks);
maxy1=log10(multpar_ks(10^-1.0978,km, cs0, cb0,Imic));
fn = @(a) (10^maxy1 - 10^a) - my*(maxy1-a)*log(10);
miny=fzero(fn, -11);

r=-10.2609 + (-8.3 - -10.269)*rand(10000,1);
y=10.^r;
mean(y)
ksm1=reshape(y, [100,100]);


fn = @(a) (10^-8.8 - 10^a) - my*(-8.8-a)*log(10);
rr2= fzero(fn, -10)

r=-9.2 + (-8.8 --9.2)*rand(10000,1);
y=10.^r;
mean(y)
ksm2=reshape(y, [100,100]);

close all
Nbin=20;
[n1, ed]=histcounts(ksm1,Nbin );
b_center1=ed(1:end-1)+diff(ed)/2;
d1=diff(ed);
figure
bar(b_center1, n1/(sum(n1)*d1(1)));   hold on

[n2, ed]=histcounts(ksm2, Nbin);
b_center2=ed(1:end-1)+diff(ed)/2;
d2=diff(ed);
bar(b_center2, n2/(sum(n2)*d2(1)));    
xlabel('k_{s,mult}')
ylabel('Normalised Rel freq')
set(gca, 'FontSize',14);


figure
area(b_center1, n1/(sum(n1)*d1(1)));
hold on 
area(b_center2, n2/(sum(n2)*d2(1)));
xlabel('k_{s,mult}')
ylabel('Normalised Rel freq')
set(gca, 'FontSize',14);

% save('C:\Users\arch9809\Box Sync\Matlab\New heterogeneous\phch_het\ksm1.mat', 'ksm1')
% save('C:\Users\arch9809\Box Sync\Matlab\New heterogeneous\phch_het\ksm2.mat', 'ksm2')