%% Steps to run the MM model
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

%% Simulate homogeneous Carbon dynamic for MM kinetics
load('field_Csp.mat');
load('field_Cb.mat');
cs0=mean2(spCs);
cb0=mean2(spCb);
co20=0;
Imic=cs0/200000;    % external C input rate fgC/h
km=((km/1000).*rho.*1e-12.*1e15).*Vol_mic;     %fgC

T=24*200; nstep=200; t=linspace(0,T,nstep); % Simulation time 200 days with daily time steps
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(t)).*Imic;    % vector of external C input at each time step 

f=@(tt,c)[-ksmm *c(1)*c(2)/(km+c(1)) + kb*c(2) + interp1(t, St,tt);...
    Y*ksmm *c(1)*c(2)/(km+c(1)) - kb*c(2); 
    (1-Y)*ksmm*c(1)*c(2)/(km+c(1))];

options = odeset('Stats','off','AbsTol',1e-6,'RelTol',1e-6);    % MATLAB ode options
[tt, Ch]=ode45(f,t,c(1,:), options); %homogeneous solution 

dec0=ksmm.*Ch(:,1).*Ch(:,2)./(km+Ch(:,1));  %homogeneous rate of decomposition
Cd=Ch.*1e-12./(rho*Vmic_cc); %Time evolution of carbon mg C at microbsite
figure;plot(tt./24,Cd(:,1)); hold on; plot(tt./24,Cd(:,2)); plot(tt./24,Cd(:,3));
xlabel('Time(day)'); ylabel('(mg C/gSoil) ');
legend('Substrate','Biomass', 'CO_2');title('Homogeneous case')
figure; plot(t./24, (1-Y).*dec0.*1e-12./(rho*Vmic_cc))
xlabel('Time(day)'); ylabel('Respirataion rate (mg C/gSoil h^{-1})');
title('Homogeneous case')


%% Biophysical heterogeneity only i.e. only substrate and microbial C are spatially heterogeneous 

%positive correlation
sav_fname=strcat('results\Ph_mm_pos');
load('field_Csp.mat');
load('field_Cb.mat');
cs=spCs;
cb=spCb;
ksmm=ones(nx,ny).*ksmm;
kb=ones(nx,ny).*kb;
km=ones(nx,ny).*km;
hetero_MM(cs,cb,sav_fname,T,t,Imic, ksmm, kb,km, Y, nx, ny, options);

%negative correlation
sav_fname=strcat('results\Ph_mm_neg');
load('field_Csn.mat');
load('field_Cb.mat');
cs=spCs;
cb=spCb;
ksmm=ones(nx,ny).*ksmm;
kb=ones(nx,ny).*kb;
km=ones(nx,ny).*km;
hetero_MM(cs,cb,sav_fname,T,t,Imic, ksmm, kb,km, Y, nx, ny, options);

%zero correlation
sav_fname=strcat('results\Ph_mm_nocorr');
load('field_Cs_noCorr.mat');
load('field_Cb.mat');
cs=spCs;
cb=spCb;
ksmm=ones(nx,ny).*ksmm;
kb=ones(nx,ny).*kb;
km=ones(nx,ny).*km;
hetero_MM(cs,cb,sav_fname,T,t,Imic, ksmm, kb,km, Y, nx, ny, options);


%% Biophysical and biochemical heterogeneity only i.e. only substrate, 
%  microbial C and kinetic parameters are spatially heterogeneous 

%positive correlation
sav_fname=strcat('results\Phch_mm_pos');
load('field_Csp.mat');
load('field_Cb.mat');
load('ksmmh1.mat');
load('kmh1.mat');
cs=spCs;
cb=spCb;
ksmm= ksmmh1;
kb=ones(nx,ny).*kb;
km=kmh1;
hetero_MM( cs,cb,sav_fname,T,t,Imic, ksmm, kb,km, Y, nx, ny, options);

%negative correlation
sav_fname=strcat('results\Phch_mm_neg');
load('field_Csn.mat');
load('field_Cb.mat');
load('ksmmh1.mat');
load('kmh1.mat');
cs=spCs;
cb=spCb;
ksmm= ksmmh1;
kb=ones(nx,ny).*kb;
km=kmh1;
hetero_MM( cs,cb,sav_fname,T,t,Imic, ksmm, kb,km, Y, nx, ny, options);


%zero correlation
sav_fname=strcat('results\Phch_mm_noCorr');
load('field_Cs_noCorr.mat');
load('field_Cb.mat');
load('ksmmh1.mat');
load('kmh1.mat');
cs=spCs;
cb=spCb;
ksmm= ksmmh1;
kb=ones(nx,ny).*kb;
km=kmh1;
hetero_MM(cs,cb,sav_fname,T,t,Imic, ksmm, kb,km, Y, nx, ny, options);
%%
rmpath(selpath)
rmpath([selpath,'\Scenario2_transientIC\Spatial_field'])
rmpath(genpath([selpath,'\Third_party_scripts\']))


