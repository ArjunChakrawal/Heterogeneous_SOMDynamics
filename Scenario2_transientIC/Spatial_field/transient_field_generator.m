% This code generates the heterogeneous fields for positive, negative and
% zero correlated substrate and microbial C for transient initial condition
% scenario. Script should be run several times until satisfactory correlations are
% found.
% IMPORTANT: For reproducebility of results, provided ksmmh1.mat, kmh1.mat and ksm1.mat
% files should be used. This script, however, includes the algorithm used
% to generate these files.

%%
clc
close all
clearvars;
selpath = uigetdir;
addpath(selpath)
addpath([selpath,'\Scenario2_transientIC\Spatial_field'])

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

faccs=12.1212;  % initial SOM fraction
total_allowedCs= soil_domain*faccs*0.01*1e15; % total amount of SOM in domain
total_allowedCb= total_allowedCs*0.01; % 1% of total substrate is microbial C


%% Generating spatially correlated Cs adn Cb with negative correlation
sp1=spatialPattern([100,100],-3);
mu=total_allowedCb/nx/ny;
sigma=10*mu;
spCb = sp1.*sigma + ones(nx,ny).*mu;
id1 =find(spCb>=0.05*max_allowedC_ms);
spCb(id1)=0.05*max_allowedC_ms;

id2 =spCb<=5000;
spCb(id2)=0;
spCb=(spCb./sum(spCb(:))).*total_allowedCb;

sp2=spatialPattern([100,100],-3);
sp3= -0.5.*spCb +  0.75*max(spCb(:)) + sp2.*max(spCb(:))*0.1;
spCs=sp3./sum(sp3(:)).*total_allowedCs ;

save('field_Csn.mat', 'spCs')
save('field_Cb.mat', 'spCb')

subplot(3,3,1)
surf(spCs);  shading flat; view(0,90); colorbar; 
title(sprintf('total Cs= %1.2d fgC',sum(spCs(:)))); axis square; axis tight;
set(gca,'FontSize',10);
subplot(3,3,2)
surf(spCb);  shading flat; view(0,90);colorbar;
title(sprintf('total Cb= %1.2d fgC',sum(spCb(:))));axis square;axis tight;
set(gca,'FontSize',10);
subplot(3,3,3)
plot(spCb(:),spCs(:),'.');xlabel('Cb');ylabel('Cs');axis square;
txt=text(max(spCb(:))/2,max(spCs(:))*0.9, sprintf('Corr(Cs,Cb)= %0.3f',corr2(spCs,spCb)));
txt.FontSize=10;txt.FontWeight='bold';
set(gca,'FontSize',10);
colormap jet

%% Generating spatially correlated Cs adn Cb with positive correlation
sp2=spatialPattern([200,200],-3);
sp2=sp2(1:100,1:100);
sp3= 200.*spCb  + max(spCb(:))*100 +(sp2).*100e6; % high Cs
sp3=sp3./sum(sp3(:)); 
spCs=sp3.*total_allowedCs ;

save('field_Csp.mat', 'spCs')

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
sp3=spatialPattern([200,200],-3);
sp3=sp3(1:100,1:100);
mu=total_allowedCs/nx/ny;
sigma=20*mu;
sp3 = sp3.*sigma + ones(nx,ny).*mu.*2;
id2 =find(sp3<=0);
sp3(id2)=0;
spCs=sp3./sum(sp3(:)).*total_allowedCs ;

save('field_Cs_noCorr.mat', 'spCs')

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


%% Generation of random field of kinetic parameters Km and ksmm for MM kinetics

r=-3 + (-1.0978 +3)*rand(10000,1);
y=10.^r;
my=mean(y);
ksmmh1=reshape(y, [nx,ny]);

kmh = 0.25 + (49.75- 0.25)*rand(10000,1);
kmh=reshape(kmh, [nx,ny]);
kmh1=(kmh.*(rho/1000).*1e-12).*1e15.*Vol_mic; %fg C

% save('ksmmh1.mat', 'ksmmh1')
% save('kmh1.mat', 'kmh1')

%% Generation of random field of kinetic parameters ks for Mult kinetics
r1=-10.1 + (-8.56 - -10.1)*rand(10000,1);
y1=10.^r1;
ksm1=reshape(y1, [100,100]);
mean(y1)
% save('ksm1.mat', 'ksm1')
%%
rmpath(selpath)
rmpath([selpath,'\Scenario2_transientIC\Spatial_field'])

