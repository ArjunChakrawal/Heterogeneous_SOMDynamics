%% Script for the Multiplicative model for scenario1 i.e. steady state initial condition

clc;close all;clearvars;
%Add path
selpath = uigetdir;
addpath(selpath)
addpath(genpath([selpath,'\Third_party_scripts\']))
addpath([selpath,'\Scenario1_SteadyStateIC\Spatial_field'])
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
Vdomain_cc=Vol_mic*nx*ny*1e-12 ; % cm^3
soil_domain=(rho*Vdomain_cc);

%% Estimation of kinetic parameter for multiplicative model ks
faccs=12.121212;  % initial SOM fraction 
total_Cs= soil_domain*faccs*0.01*1e15; % total amount of SOM in domain
total_Cb= total_Cs*0.01; % 1% of total substrate is microbial C
cs0=(total_Cs/(nx*ny))*1e-12./(rho*Vmic_cc);cb0=(total_Cb/(nx*ny))*1e-12./(rho*Vmic_cc);co20=0;
T=24*1000;Imic=cs0/200000;
ks=ks_mult_estimation(T,cs0,cb0,co20,kb,ksmm, km, Imic,Y,rho,Vol_mic);
%% Steady state substrate and microbial C 
Imic=(total_Cs/(nx*ny))/200000;
css=kb/(Y*ks);      %Steady state substrate
cbs=Y*Imic/(kb*(1-Y));      %Steady state microbial C 

%% Homogeneous case SOM dynamics
cs0=css;cb0=cbs;co20=0;
options = odeset('Stats','off','AbsTol',1e-6,'RelTol',1e-6);
T=24*20*365; nstep=200; t=linspace(0,T,nstep);
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(t)).*Imic;

[t1,C1]= DEC_multiplicative(cs0,cb0,co20,t, ks, kb, Y,St, options);
dec0=ks.*C1(:,1).*C1(:,2);
C1=C1.*1e-12./(rho*Vmic_cc); %Time evolution of carbon mg C/gSoil

figure;
plot(t1./24./365,C1(:,1)); hold on; plot(t1./24./365,C1(:,2));
xlabel('Time(year)'); ylabel('(mg C/gSoil) ');
legend('Substrate','Biomass');title('Homogeneous case')

%% Biophysical heterogeneity only i.e. only substrate and microbial C are spatially heterogeneous 

%positive correlation
sav_fname=strcat('results\Ph_mult_pos_ss_1');
load('ph_pos_cs.mat');
load('ph_cb.mat');
cs=spCs;
cb=spCb;
ks=ones(nx,ny).*ks;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%negative correlation
sav_fname=strcat('results\Ph_mult_neg_ss_1');
load('ph_neg_cs.mat');
load('ph_cb.mat');
cs=spCs;
cb=spCb;
ks=ones(nx,ny).*ks;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%zero correlation
sav_fname=strcat('results\Ph_mult_zero_ss_1');
load('ph_zero_cs.mat');
load('ph_cb.mat');
cs=spCs;
cb=spCb;
ks=ones(nx,ny).*ks;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%% Biophysical and biochemical heterogeneity only i.e. only substrate, 
%  microbial C and kinetic parameters are spatially heterogeneous 

%positive correlation
sav_fname=strcat('results\Ph1ch2_mult_pos_ss');
load('ph_pos_cs.mat');
load('ph_cb.mat');
load('ksm2.mat');
cs=spCs;
cb=spCb;
ks=ksm2;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%negative correlation
sav_fname=strcat('results\Ph1ch2_mult_neg_ss');
load('ph_neg_cs.mat');
load('ph_cb.mat');
load('ksm2.mat');
cs=spCs;
cb=spCb;
ks=ksm2;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic, ks, kb, Y, nx, ny, options)

%zero correlation
sav_fname=strcat('results\Ph1ch2_mult_zero_ss');
load('ph_zero_cs.mat');
load('ph_cb.mat');
load('ksm2.mat');
cs=spCs;
cb=spCb;
ks=ksm2;
kb=ones(nx,ny).*kb;
hetero_Mult(cs,cb,sav_fname,T,t,Imic,ks, kb, Y, nx, ny, options)

%%
rmpath(selpath)
rmpath([selpath,'\Scenario1_SteadyStateIC\Spatial_field'])
rmpath(genpath([selpath,'\Third_party_scripts\']))

