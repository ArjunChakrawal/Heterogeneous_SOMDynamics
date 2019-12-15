clc;
close all;
clearvars;

selpath = uigetdir;
addpath(selpath)
addpath(genpath([selpath,'\Third_party_scripts\']))
addpath([selpath,'\Scenario2_transientIC\Mult_transport\'])


cases = 2;
sig = 1:0.5:10;
ix = [1, 4, 9, 11, 15, 19];
Nl = length(ix) + 1;
LC = linspecer(Nl); % for disctinctive line color

ksmm = 0.018;
kb = 0.00028; %h-1
Y = 0.31;
km = 25;
rho = 1.65; %g/cm3 soil bulk density
rho_wood = 1.1; %g/cm3  max carbon density in form of wood
nx = 100;
ny = 100; %One pixel = 25micron

Vol_mic = 50 * 50 * 50; %um^3
Vmic_cc = Vol_mic * 1e-12; %cm3
Vdomain_cc = Vol_mic * nx * ny * 1e-12; % cm^3
T = 24 * 200; % Simulation time 200 days with daily time steps
dt = 24;
Nt = round(T/dt) + 1;
ttime = 0:dt:T; % Mesh points in time

%% homo
load('output\Ph_mult_pos_0.mat')
cs0=mean2(spCs);cb0=mean2(spCb);co20=0;
Imic=cs0/200000;
ks = mean2(ks);
kb = mean2(kb);
options = odeset('Stats','off','AbsTol',1e-6,'RelTol',1e-6);
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(ttime)).*Imic;

c(1,:)=[cs0,cb0,co20];
f=@(tt,c)[-ks *c(1)*c(2) + kb*c(2) + interp1(ttime, St,tt);...
    Y*ks *c(1)*c(2) - kb*c(2); (1-Y)*ks*c(1)*c(2)];
[~, Chh]=ode45(f,ttime,c(1,:), options);
dec0h = ks .* Chh(:, 1) .* Chh(:, 2) ; %homogeneous rate of decomposition


%% Figure properties
phi = (1+sqrt(5))/2; %godlen ratio
width =7;     % Width in inches
height = 8;    % Height in inches
alw = 0.5;    % AxesLineWidth
fsz = 12;      % FontSize
fszax=fsz-1;
fsz_title= fsz+1;
fsz_ax_lgnd = fsz;
lws=1.75;
lw = 0.75;      % LineWidth
lwdot= 0.5;
lwdotmfa= 1.75;
xmsz=5;
msz =3.5;       % MarkerSize
clr = [0.1,0.1,0.1];     % Line color
xclr=[0.3,0.3,0.3];
ftype1= 'Helvetica';  %hLegend, gca
ftype2='Times New Roman'; % hTitle, hXLabel, hYLabel
xycolor2= [0.1 0.1 0.1]*0;
fig=figure;
pos = get(gcf, 'Position');
set(gcf, 'Position', [20 10 width*100, height*100]); %<- Set size
% [ha, pos] = tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
[ha, pos] = tight_subplot(3,2, [.05 .05],[.225 .05],[.2 .05]);
ax1=ha(1);a2=ha(3);a3=ha(5);
a4=ha(2);a5=ha(4);a6=ha(6);

TT=60;
xtick= 0:30:TT;
cstck=[0,40,80];
%% Post-processing Positive
load('output\Ph_mult_pos_0.mat')
plot(ax1,ttime./24, (1 - Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws);hold(ax1,'on'); 
plot(ax1,ttime./24,(1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa);
plot(ax1,ttime./24, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws );  
plot(ax1,ttime./24, (1-Y).*dec0h.*1e-12./(rho*Vmic_cc),'-','linewidth',lws);  


load('output\Ph_mult_neg_0.mat')
plot(a2,ttime./24, (1 - Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws);hold(a2,'on'); 
plot(a2,ttime./24,(1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa);
plot(a2,ttime./24, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws );  
plot(a2,ttime./24, (1-Y).*dec0h.*1e-12./(rho*Vmic_cc),'-','linewidth',lws);  


load('output\Ph_mult_noCorr_0.mat')
plot(a3,ttime./24, (1 - Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws);hold(a3,'on'); 
plot(a3,ttime./24,(1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa);
plot(a3,ttime./24, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws );  
plot(a3,ttime./24, (1-Y).*dec0h.*1e-12./(rho*Vmic_cc),'-','linewidth',lws);  

%%
tmp=1:3:200;

load('output\Ph_mult_fullHetero_pos_0.mat')
cov1=(1-Y).*cov_csks.*Cb_mean;
cov2=(1-Y).*cov_cbks.*Cs_mean;
cov3=(1-Y).*cov_cscb;
cov4=(1-Y).*E_cs_cb_ks;

p1=plot(a4,ttime./24, (1 - Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws);hold(a4,'on');
p2=plot(a4,ttime./24,(1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa );
p3=plot(a4,ttime./24, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws );  
p4=plot(a4,ttime./24, (1-Y).*dec0h.*1e-12./(rho*Vmic_cc),'-','linewidth',lws);  


p5a=plot(a4,t(tmp)./24,cov1(tmp).*1e-12./(rho*Vmic_cc),'-x','linewidth',lwdot,'MarkerSize',xmsz,'MarkerEdgeColor',xycolor2); 
p6a=plot(a4,t(tmp)./24,cov2(tmp).*1e-12./(rho*Vmic_cc),'-o','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor2); 
p7a=plot(a4,t(tmp)./24,cov3(tmp).*mean2(ks).*1e-12./(rho*Vmic_cc),'-s','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor2); 
p8a=plot(a4,t(tmp)./24,cov4(tmp).*1e-12./(rho*Vmic_cc),'-d','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor2); 

load('output\Ph_mult_fullHetero_neg_0.mat')

cov1=(1-Y).*cov_csks.*Cb_mean;
cov2=(1-Y).*cov_cbks.*Cs_mean;
cov3=(1-Y).*cov_cscb;
cov4=(1-Y).*E_cs_cb_ks;
p1=plot(a5,ttime./24, (1 - Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws);hold(a5,'on'); 
p2=plot(a5,ttime./24,(1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa );
p3=plot(a5,ttime./24, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws );  
p4=plot(a5,ttime./24, (1-Y).*dec0h.*1e-12./(rho*Vmic_cc),'-','linewidth',lws);  

p5a=plot(a5,t(tmp)./24,cov1(tmp).*1e-12./(rho*Vmic_cc),'-x','linewidth',lwdot,'MarkerSize',xmsz,'MarkerEdgeColor',xycolor2); 
p6a=plot(a5,t(tmp)./24,cov2(tmp).*1e-12./(rho*Vmic_cc),'-o','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor2); 
p7a=plot(a5,t(tmp)./24,cov3(tmp).*mean2(ks).*1e-12./(rho*Vmic_cc),'-s','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor2); 
p8a=plot(a5,t(tmp)./24,cov4(tmp).*1e-12./(rho*Vmic_cc),'-d','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor2); 


load('output\Ph_mult_fullHetero_noCorr_0.mat')
cov1=(1-Y).*cov_csks.*Cb_mean;
cov2=(1-Y).*cov_cbks.*Cs_mean;
cov3=(1-Y).*cov_cscb;
cov4=(1-Y).*E_cs_cb_ks;
p1=plot(a6,ttime./24, (1 - Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws);hold(a6,'on'); 
p2=plot(a6,ttime./24,(1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa );
p3=plot(a6,ttime./24, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws );  
p4=plot(a6,ttime./24, (1-Y).*dec0h.*1e-12./(rho*Vmic_cc),'-','linewidth',lws);  

p5a=plot(a6,t(tmp)./24,cov1(tmp).*1e-12./(rho*Vmic_cc),'-x','linewidth',lwdot,'MarkerSize',xmsz,'MarkerEdgeColor',xycolor2); 
p6a=plot(a6,t(tmp)./24,cov2(tmp).*1e-12./(rho*Vmic_cc),'-o','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor2); 
p7a=plot(a6,t(tmp)./24,cov3(tmp).*mean2(ks).*1e-12./(rho*Vmic_cc),'-s','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor2); 
p8a=plot(a6,t(tmp)./24,cov4(tmp).*1e-12./(rho*Vmic_cc),'-d','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor2); 

%%
set(ax1,'XLim',[0,TT],'XTickLabel',[],'YLim',[-0.15,0.2],'YTick',-0.2:0.1:0.2,...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor2, 'YColor', xycolor2);

set(a2,'XLim',[0,TT],'XTickLabel',[],'YLim',[-0.15,0.2],'YTick',-0.2:0.1:0.2,...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off', ...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor2, 'YColor', xycolor2);

set(a3,'XLim',[0,TT],'YLim',[-0.15,0.2],'YTick',-0.2:0.1:0.2,...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off', ...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor2, 'YColor', xycolor2);

set(a4,'XLim',[0,TT],'XTickLabel',[],'YTickLabel',[],'YLim',[-0.15,0.2],'YTick',-0.2:0.1:0.2,...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor2, 'YColor', xycolor2);
% axis(ha(1), 'square');
set(a5,'XLim',[0,TT],'XTickLabel',[],'YTickLabel',[],'YLim',[-0.15,0.2],'YTick',-0.2:0.1:0.2,...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor2, 'YColor', xycolor2);
% axis(ha(2), 'square');
set(a6,'XLim',[0,TT],'YTickLabel',[],'YLim',[-0.15,0.2],'YTick',-0.2:0.1:0.2,...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor2, 'YColor', xycolor2);

%%
hxlable=xlabel(a6,'Time (day)','FontSize', fsz, 'FontName', ftype2,'Interpreter','latex');
hx=hxlable.Position;
hxlable.Position=[-5.0000   -0.2146   -1.0000];
% [hx(1)+23,hx(2)-0.00001,hx(3)];

ylabel(ax1,'$\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a2,'$\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a3,'$\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
align_Ylabels(gcf)

axes(a6); %set the current axes to axes2
[cl,hobj]=columnlegend(2,{'$\overline{R}_{het}$', 'MFA','$ \sum $ HOT','$\overline{R}_{hom}$',...
    '$(1-Y)\overline{C}_b$ $\overline{k_{s,mult}^{\prime} C_s^{\prime}}$',...
    '$(1-Y)\overline{C}_s$ $\overline{k_{s,mult}^{\prime} C_b^{\prime}}$',...
    '$(1-Y)\overline{k}_{s,mult}$ $\overline{C_s^{\prime} C_b^{\prime}}$',...
    '$(1-Y)\overline{k_{s,mult}^{\prime} C_s^{\prime} C_b^{\prime}}$'}, 'location','northeast','Interpreter','latex');
objhl = findobj(hobj, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 6); %// set marker size as desired
cl.Position=[0.3330   -0.0371    0.4354    0.1981];

axes(ax1);
t1=text(-20,-.05,'Positive','Interpreter','latex','FontSize', fsz, 'FontName', ftype2, 'FontWeight',   'bold' );
t1.Rotation=90;

axes(a2);
t1=text(-20,-.05,'Negative','Interpreter','latex','FontSize', fsz, 'FontName', ftype2, 'FontWeight',   'bold' );
t1.Rotation=90;
axes(a3);
t1=text(-20,-.05,'Uncorrelated','Interpreter','latex','FontSize', fsz, 'FontName', ftype2, 'FontWeight',   'bold' );
t1.Rotation=90;

ta1=title(ax1, 'Biophysical heterogeneity','FontSize', fsz_title, 'FontName', ftype2,'FontWeight' , 'normal' );
tp1=get(ta1, 'Position');
set(ta1, 'Position', [tp1(1), 0.225])
ta2=title(a4, 'Full heterogeneity','FontSize', fsz_title, 'FontName', ftype2,'FontWeight' , 'normal');
tp2=get(ta2, 'Position');
set(ta2, 'Position', [tp2(1),0.225])
str={'a','b','c'};k=1;
for i=1:2:6
    axes(ha(i)); %set the current axes to axes2
    text(TT-5, ha(i).YLim(2) - ha(i).YLim(2)/10, ['(',str{k},')'], ...
        'FontWeight','normal','FontSize', fsz, 'FontName', ftype2)
    k=k+1;
end
str={'d','e','f'};
k=1;
for i=2:2:6
    
    axes(ha(i)); %set the current axes to axes2
    text(TT-5, ha(i).YLim(2) - ha(i).YLim(2)/10, ['(',str{k},')'], ...
        'FontWeight','normal','FontSize', fsz, 'FontName', ftype2)
    k=k+1;
end

set(gcf, 'Color','w');
export_fig(gcf, 'FigA4.pdf', '-r300')

%%
rmpath(selpath)
rmpath(genpath([selpath,'\Third_party_scripts\']))
rmpath(genpath([selpath,'\Scenario2_transientIC\Mult_transport\']))