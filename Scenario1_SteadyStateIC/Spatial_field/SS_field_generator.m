% This code generates the heterogeneous fields for positive, negative and
% zero correlated substrate and microbial C for steadt state initial condition
% scenario. Script should be run several times until satisfactory correlations are
% found.
% IMPORTANT: For reproducebility of results, provided ksm1.mat
% file should be used. This script, however, includes the algorithm used
% to generate these files.
%%
clc
close all
clearvars;
selpath = uigetdir;
addpath(selpath)
addpath(genpath([selpath,'\Third_party_scripts\']))
addpath([selpath,'\Scenario1_SteadyStateIC\Spatial_field\'])

ksmm=0.018;
kb=0.00028;  %h-1
Y=0.31;
km=25;
rho= 1.65;  %g/cm3 soil bulk density 
rho_SOM=1.1;    %g/cm3  SOM density
nx=100;ny=100;  %One pixel = 25micron

Vol_mic = 50*50*50; %um^3
Vmic_cc=Vol_mic*1e-12; %cm3
Vdomain_cc=Vol_mic*nx*ny*1e-12 ; % cm^3
soil_domain=(rho*Vdomain_cc);

maxC_ms=(rho_SOM*1e-12)*Vol_mic*1e15;   %fg C
max_allowedC_ms=maxC_ms.*0.5; % approx. 50% of SOM is organic C

%% Estimation of kinetic parameter for multiplicative model ks
faccs=12.121212;  % initial SOM fraction 
total_Cs= soil_domain*faccs*0.01*1e15; % total amount of SOM in domain
total_Cb= total_Cs*0.01; % 1% of total substrate is microbial C
cs0=(total_Cs/(nx*ny))*1e-12./(rho*Vmic_cc);cb0=(total_Cb/(nx*ny))*1e-12./(rho*Vmic_cc);co20=0;
T=24*1000;Imic=cs0/200000;
ks0=ks_mult_estimation(T,cs0,cb0,co20,kb,ksmm, km, Imic,Y,rho,Vol_mic);

Imic=(total_Cs/(nx*ny))/200000;

%% Steady state substrate and microbial C 
ks=ks0;
css=kb/(Y*ks);      %Steady state substrate
cbs=Y*Imic/(kb*(1-Y));      %Steady state microbial C
total_allowedCs= css*(nx*ny);   % total amount of SOM in domain
total_allowedCb= cbs*(nx*ny);   % total amount of microbial C  in domain
load('ph_cb.mat');

sp2=spatialPattern([100,100],-3);
sp3= 200.*spCb  + max(spCb(:))*100 +(sp2).*5e7; % high Cs
sp3=sp3./sum(sp3(:)); 
spCs=sp3./sum(sp3(:)).*total_allowedCs ;

% save('ph_pos_cs_ks1.mat', 'spCs')


figure
subplot(3,3,1)
surf(spCs);  shading flat; view(0,90); colorbar; title(sprintf('total Cs= %1.2d fgC',sum(spCs(:)))); 
axis square; axis tight;
set(gca,'FontSize',10);
subplot(3,3,2)
surf(spCb);  shading flat; view(0,90);colorbar;title(sprintf('total Cb= %1.2d fgC',sum(spCb(:))));
axis square;axis tight;
set(gca,'FontSize',10);
subplot(3,3,3)
plot(spCb(:),spCs(:),'.');xlabel('Cb');ylabel('Cs');axis square;
txt=text(max(spCb(:))/2,max(spCs(:))*0.9, sprintf('Corr(Cs,Cb)= %0.3f',corr2(spCs,spCb)));
txt.FontSize=12;txt.FontWeight='bold';
set(gca,'FontSize',10);
colormap jet
%% Generating spatially correlated Cs adn Cb with positive correlation
ks=1/2*ks0;
css=kb/(Y*ks);      %Steady state substrate
total_allowedCs= css*(nx*ny);   % total amount of SOM in domain

sp2=spatialPattern([100,100],-3);
sp3= 200.*spCb  + max(spCb(:))*100 +(sp2).*5e7; % high Cs
sp3=sp3./sum(sp3(:)); 
spCs=sp3./sum(sp3(:)).*total_allowedCs ;

% save('ph_pos_cs_ks2.mat', 'spCs')

subplot(3,3,4)
surf(spCs);  shading flat; view(0,90); colorbar; 
title(sprintf('total Cs= %1.2d fgC',sum(spCs(:)))); axis square; axis tight;
set(gca,'FontSize',10);
subplot(3,3,5)
surf(spCb);  shading flat; view(0,90);colorbar;
title(sprintf('total Cb= %1.2d fgC',sum(spCb(:))));axis square;axis tight;
set(gca,'FontSize',10);
subplot(3,3,6)
plot(spCb(:),spCs(:),'.');xlabel('Cb');ylabel('Cs');axis square;
txt=text(max(spCb(:))/2,max(spCs(:))*0.9, sprintf('Corr(Cs,Cb)= %0.3f',corr2(spCs,spCb)));
txt.FontSize=10;txt.FontWeight='bold';
set(gca,'FontSize',10);
%% Generating spatially correlated Cs adn Cb with zero correlation
ks=1/4 *ks0;
css=kb/(Y*ks);      %Steady state substrate
total_allowedCs= css*(nx*ny);   % total amount of SOM in domain

sp2=spatialPattern([100,100],-3);
sp3= 200.*spCb  + max(spCb(:))*100 +(sp2).*5e7; % high Cs
sp3=sp3./sum(sp3(:)); 
spCs=sp3./sum(sp3(:)).*total_allowedCs ;

% save('ph_pos_cs_ks3.mat', 'spCs')

subplot(3,3,7)
surf(spCs);  shading flat; view(0,90); colorbar; 
title(sprintf('total Cs= %1.2d fgC',sum(spCs(:)))); axis square; axis tight;
set(gca,'FontSize',10);
subplot(3,3,8)
surf(spCb);  shading flat; view(0,90);colorbar;
title(sprintf('total Cb= %1.2d fgC',sum(spCb(:))));axis square;axis tight;
set(gca,'FontSize',10);
subplot(3,3,9)
plot(spCb(:),spCs(:),'.');xlabel('Cb');ylabel('Cs');axis square;
txt=text(max(spCb(:))/2,max(spCs(:))*0.9, sprintf('Corr(Cs,Cb)= %0.3f',corr2(spCs,spCb)));
txt.FontSize=10;txt.FontWeight='bold';
set(gca,'FontSize',10);

%% Generation of random field of kinetic parameters ks for Mult kinetics
r1=-10.1 + (-8.56 - -10.1)*rand(10000,1);
y1=10.^r1;
ksm1=reshape(y1, [100,100]);
mean(y1)
% save('ksm1.mat', 'ksm1')
r2=-9.4 + (-8.9 --9.4)*rand(10000,1);
y2=10.^r2;
mean(y2)
ksm2=reshape(y2, [100,100]);
% save('ksm2.mat', 'ksm2')

% a= -10.1;b=-8.56;
% (0.00028 *(10^-a - 10^-b)/(0.31*(b-a)*log(10))) .*1e-12./(rho*Vmic_cc)
% 
% a= -9.4;b=-8.9;
% (0.00028 *(10^-a - 10^-b)/(0.31*(b-a)*log(10)) ) .*1e-12./(rho*Vmic_cc)

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
%%
rmpath(selpath)
rmpath([selpath,'\Scenario1_SteadyStateIC\Spatial_field\'])
rmpath(genpath([selpath,'\Third_party_scripts\']))
