clc
close all
clearvars;
warning('off', 'all')

selpath = uigetdir;
addpath(selpath)
addpath(genpath([selpath,'\Third_party_scripts\']))
addpath([selpath,'\Scenario2_transientIC\'])

cases=2;
sig=1:0.5:10; 
ix=[1,4,9,11,15,19];
Nl=length(ix)+1;
LC=linspecer(Nl); % for disctinctive line color 

ksmm=0.018;
kb=0.00028;  %h-1
Y=0.31;
km=25;
rho= 1.65; %g/cm3 soil bulk density 
nx=100;ny=100; %One pixel = 25micron

Vol_mic = 50*50*50; %um^3
Vmic_cc=Vol_mic*1e-12; %cm3
%% Figure properties
phi = (1+sqrt(5))/2; %godlen ratio
width = 5;     % Width in inches
height = 5.25;    % Height in inches
alw = 0.5;    % AxesLineWidth
fsz = 10;      % FontSize
fsz_title= fsz+1;
fsz_ax_lgnd = fsz-1;
lw = 1.0;      % LineWidth
msz = 8;       % MarkerSize
clr = 'k';     % Line color
ftype1= 'Helvetica';  %hLegend, gca
ftype2='Times New Roman'; % hTitle, hXLabel, hYLabel
xycolor= [0.2 0.2 0.2];
fig=figure;
pos = get(gcf, 'Position');
set(gcf, 'Position', [1228 267 width*100, height*100]); %<- Set size
% [ha, pos] = tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
[ha, pos] = tight_subplot(3, 2, [.035 .05],[.25 .1],[.15 .05]);
ax1=ha(1);a2=ha(2);a3=ha(3);a4=ha(4);
a5=ha(5); a6=ha(6);
%% homo mult
load('Mult_transport\output\Ph_mult_pos_0.mat')
cs0=mean2(spCs);cb0=mean2(spCb);co20=0;
ks=mean2(ks);
kb=mean2(kb);
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(t)).*Imic;
[~,C1]= DEC_multiplicative(cs0,cb0,co20,t, ks, kb, Y,St, options);
dec_homo=ks.*C1(:,1).*C1(:,2);


load('Mult_transport\output\Ph_mult_pos_0.mat')
plot(ax1,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'--','linewidth',lw, 'Color', clr);  hold(ax1,'on');
load('Mult_transport\output\Ph_mult_neg_0.mat')
plot(ax1,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'-.','linewidth',lw, 'Color', clr);   
load('Mult_transport\output\Ph_mult_noCorr_0.mat')
plot(ax1,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),':','linewidth',lw, 'Color', clr); 
plot(ax1,C1(:,1).*1e-12./(rho*Vmic_cc), ((1-Y).*dec_homo./C1(:,2)),'-','linewidth',lw, 'Color', clr); 

load('Mult_transport\output\Ph_mult_fullHetero_pos_0.mat')
plot(a2,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'--','linewidth',lw, 'Color', clr);  hold(a2,'on');
load('Mult_transport\output\Ph_mult_fullHetero_neg_0.mat')
plot(a2,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'-.','linewidth',lw, 'Color', clr);   
load('Mult_transport\output\Ph_mult_fullHetero_noCorr_0.mat')
plot(a2,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),':','linewidth',lw, 'Color', clr);  
plot(a2,C1(:,1).*1e-12./(rho*Vmic_cc), ((1-Y).*dec_homo./C1(:,2)),'-','linewidth',lw, 'Color', clr); 

%% 
load('MM_transient\results\output_transport\varing_alpha\results_v3_Ph_mm_pos_0.mat')
cs0=mean2(spCs);cb0=mean2(spCb);co20=0;
ksmm=mean2(ksmm);
km=mean2(km);
kb=mean2(kb);
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(t)).*Imic;
f=@(tt,c)[-ksmm *c(1)*c(2)/(km+c(1)) + kb*c(2) + interp1(t, St,tt);...
    Y*ksmm *c(1)*c(2)/(km+c(1)) - kb*c(2); (1-Y)*ksmm*c(1)*c(2)/(km+c(1))];
[~, C1]=ode45(f,t,c(1,:), options); %homogeneous solution 
dec_homo=ksmm.*C1(:,1).*C1(:,2)./(km+C1(:,1));


load('MM_transient\results\output_transport\varing_alpha\results_v3_Ph_mm_pos_0.mat');
plot(a3,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'--','linewidth',lw, 'Color', clr);  hold(a3,'on');
load('MM_transient\results\output_transport\varing_alpha\results_v3_Ph_mm_neg_0.mat');
plot(a3,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'-.','linewidth',lw, 'Color', clr);   
load('MM_transient\results\output_transport\varing_alpha\results_v3_Ph_mm_noCorr_0.mat');
plot(a3,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),':','linewidth',lw, 'Color', clr);  
plot(a3,C1(:,1).*1e-12./(rho*Vmic_cc), ((1-Y).*dec_homo./C1(:,2)),'-','linewidth',lw, 'Color', clr); 

load('MM_transient\results\output_transport\varing_alpha\results_v3_Ph_mm_fullHetero_pos_0.mat');
plot(a4,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'--','linewidth',lw, 'Color', clr);  hold(a4,'on');
load('MM_transient\results\output_transport\varing_alpha\results_v3_Ph_mm_fullHetero_neg_0.mat');
plot(a4,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'-.','linewidth',lw, 'Color', clr);   
load('MM_transient\results\output_transport\varing_alpha\results_v3_Ph_mm_fullHetero_noCorr_0.mat');
plot(a4,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),':','linewidth',lw, 'Color', clr);  
plot(a4,C1(:,1).*1e-12./(rho*Vmic_cc), ((1-Y).*dec_homo./C1(:,2)),'-','linewidth',lw, 'Color', clr); 


%%
mm_in_path='MM_Inv_transient\output\';
load([mm_in_path,'Ph_mm_inv_pos_0.mat'])
plot(a5,Cs_mean.*1e-12./(rho*Vmic_cc), (1-Y).*dec_m./Cb_mean,'--','linewidth',lw, 'Color', clr);  hold(a5,'on');
load([mm_in_path,'Ph_mm_inv_neg_0.mat'])
plot(a5,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'-.','linewidth',lw, 'Color', clr);   
load([mm_in_path,'Ph_mm_inv_noCorr_0.mat'])
plot(a5,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),':','linewidth',lw, 'Color', clr);  
plot(a5,Ch2(:,1).*1e-12./(rho*Vmic_cc), (1-Y).*decmminv./Ch2(:,2),'-','linewidth',lw, 'Color', clr); 

load([mm_in_path,'Ph_mm_inv_fullhetero_pos_0.mat'])
plot(a6,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'--','linewidth',lw, 'Color', clr);  hold(a6,'on');
load([mm_in_path,'Ph_mm_inv_fullhetero_neg_0.mat'])
plot(a6,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),'-.','linewidth',lw, 'Color', clr);   
load([mm_in_path,'Ph_mm_inv_fullhetero_noCorr_0.mat'])
plot(a6,Cs_mean.*1e-12./(rho*Vmic_cc), ((1-Y).*dec_m./Cb_mean),':','linewidth',lw, 'Color', clr);  
plot(a6,Ch2(:,1).*1e-12./(rho*Vmic_cc), (1-Y).*decmminv./Ch2(:,2),'-','linewidth',lw, 'Color', clr); 

%%
hxlable=xlabel(a6,'$\overline{C}_S$ (mgC/gSoil)','Interpreter','latex','FontSize', fsz, 'FontName', ftype2);
hx=hxlable.Position;
hxlable.Position= [hx(1)-90,hx(2)-0.006,hx(3)];

hylabel1=ylabel(a3,'$SGR$ \ $\mathrm{(h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2);
% hy=hylabel1.Position;
% hylabel1.Position=[hy(1),hy(2)-0.018,hy(3)];
hlegend=legend(ax1,'Positive', 'Negative','Uncorrelated','Homogeneous');
set(hlegend,'box', 'off','Location', 'northwest', 'Interpreter','latex','FontSize', fsz, 'FontName', ftype1);
hlegend.NumColumns=2;
hlegend.Position = [0.3151    0.0790     0.4328    0.0666];
ylim(ha,[0, 0.04]);
set(ha(1:4),'XTickLabel',[]);
set(ha(2:2:6),'YTickLabel',[]);
set(ha,'FontSize', fsz,'FontName',ftype1, 'LineWidth',  alw,...
    'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on',...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);

% Biophysical heterogeneity \newline
% Full heterogeneity
text(ax1,20,0.035, 'Multiplicative','FontSize', fsz, 'FontName', ftype2 )
text(a3,20,0.035, 'MM','FontSize', fsz, 'FontName', ftype2 )
text(a5,20,0.035, 'Inverse MM','FontSize', fsz, 'FontName', ftype2 )


ta1 = title(ax1, 'Biophysical heterogeneity', 'FontSize', fsz_title, 'FontName', ftype2, 'FontWeight', 'normal');
ta2 = title(a2, 'Full heterogeneity', 'FontSize', fsz_title, 'FontName', ftype2, 'FontWeight', 'normal');
ta1.Position = [75, 0.042, 0];
ta2.Position = ta1.Position;

numbers = 1:26;
letters = lower(char(numbers+64));
for i = 1:length(ha)
    axes(ha(i)); %set the current axes to axes2
    text(ha(i).XLim(2)-15, ha(i).YLim(2)-ha(i).YLim(2)/10, ['(', letters(i), ')'], 'FontWeight', 'normal', 'FontSize', fsz+2, 'FontName', ftype2)
end
align_Ylabels(gcf)

set(gcf, 'Color', 'w')
export_fig(gcf,'Fig8.pdf','-r300');
%%
rmpath(selpath)
rmpath(genpath([selpath,'\Third_party_scripts\']))
rmpath([selpath,'\Scenario2_transientIC\'])
warning('on', 'all')
