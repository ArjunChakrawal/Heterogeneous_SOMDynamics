%% Steps to run the Multiplicative model
% Step 1:  Run field_trans.m to generate heterogeneous fields of substrate
% and microbial C.
% Step 2: 

%%
clc;close all;clearvars;
selpath = uigetdir;
addpath(selpath)
addpath(genpath([selpath,'\Third_party_scripts\']))
addpath([selpath,'\Scenario2_transientIC\Spatial_field'])

%% define model parameters
ksmm=0.018;     % h-1
kb=0.00028;     % h-1
Y=0.31;         % microbial carbon use efficiency
km=25;          % mgC/gSoil
rho= 1.65;      % g/cm3 soil bulk density 
rho_wood=1.1;   % g/cm3  max carbon density in form of wood 
nx=100;ny=100;  % x y grid points
Vol_mic = 50*50*50;     % volume of a microsite in µm^3 (each grid cell)
Vmic_cc=Vol_mic*1e-12;  %cm3

%% ksmult estimation
load('field_Csp.mat');
load('field_Cb.mat');
cs0=mean2(spCs)*1e-12./(rho*Vmic_cc);cb0=mean2(spCb)*1e-12./(rho*Vmic_cc);co20=0;
Imic=cs0/200000; T=24*1000;
ks=ks_mult_estimation(T,cs0,cb0,co20,kb,ksmm, km, Imic,Y,rho,Vol_mic);

%% Homogeneous  case 
cs0=mean2(spCs);cb0=mean2(spCb);co20=0;
Imic=cs0/200000;
options = odeset('Stats','off','AbsTol',1e-6,'RelTol',1e-6);
T=24*200; nstep=200; t=linspace(0,T,nstep);
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(t)).*Imic;

c(1,:)=[cs0,cb0,co20];
f=@(tt,c)[-ks *c(1)*c(2) + kb*c(2) + interp1(t, St,tt);...
    Y*ks *c(1)*c(2) - kb*c(2); (1-Y)*ks*c(1)*c(2)];
[tt, Ch]=ode45(f,t,c(1,:), options);

dec0=ks.*Ch(:,1).*Ch(:,2);
Cd=Ch.*1e-12./(rho*Vmic_cc); %Time evolution of carbon mg C at microbsite
figure;plot(tt./24,Cd(:,1)); hold on; plot(tt./24,Cd(:,2)); plot(tt./24,Cd(:,3));
xlabel('Time(day)'); ylabel('(mg C/gSoil) ');
legend('Substrate','Biomass', 'CO_2');title('Homogeneous case')
figure; plot(t./24, (1-Y).*dec0.*1e-12./(rho*Vmic_cc))
xlabel('Time(day)'); ylabel('Respirataion rate (mg C/gSoil h^{-1})');
title('Homogeneous case')
%% Biophysical heterogeneity only i.e. only substrate and microbial C are spatially heterogeneous 

%positive correlation
sav_fname=strcat('results\Ph_mult_pos');
load('field_Csp.mat');
load('field_Cb.mat');
cs=spCs;
cb=spCb;
ks=ones(nx,ny).*ks;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%negative correlation
sav_fname=strcat('results\Ph_mult_neg');
load('field_Csn.mat');
load('field_Cb.mat');
cs=spCs;
cb=spCb;
ks=ones(nx,ny).*ks;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%zero correlation
sav_fname=strcat('results\Ph_mult_noCorr');
load('field_Cs_noCorr.mat');
load('field_Cb.mat');
cs=spCs;
cb=spCb;
ks=ones(nx,ny).*ks;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%% Biophysical and biochemical heterogeneity only i.e. only substrate, 
%  microbial C and kinetic parameters are spatially heterogeneous 

%positive correlation
sav_fname=strcat('results\Phch_mult_pos');
load('field_Csp.mat');
load('field_Cb.mat');
load('ksm1.mat');
cs=spCs;
cb=spCb;
ks=ksm1;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%negative correlation
sav_fname=strcat('results\Phch_mult_neg');
load('field_Csn.mat');
load('field_Cb.mat');
load('ksm1.mat');
cs=spCs;
cb=spCb;
ks=ksm1;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%zero correlation
sav_fname=strcat('results\Phch_mult_noCorr');
load('field_Cs_noCorr.mat');
load('field_Cb.mat');
load('ksm1.mat');
cs=spCs;
cb=spCb;
ks=ksm1;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic,ks, kb, Y, nx, ny, options)

%%

rmpath(selpath)
rmpath([selpath,'\Scenario2_transientIC\Spatial_field'])
rmpath(genpath([selpath,'\Third_party_scripts\']))

