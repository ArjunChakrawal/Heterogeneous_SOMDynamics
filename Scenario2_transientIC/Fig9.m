clc
close all
clearvars;

selpath = uigetdir;
addpath(selpath)
addpath(genpath([selpath,'\Third_party_scripts\']))
addpath([selpath,'\Scenario2_transientIC\'])

LC=linspecer(10); % for disctinctive line color 

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
set(gcf, 'Position', [1228 267 width*100, height*100]); %<- Set size
% [ha, pos] = tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
[ha, pos] = tight_subplot(3, 2, [.035 .05],[.25 .1],[.15 .05]);
ax1=ha(1);a2=ha(2);a3=ha(3);a4=ha(4);
a5=ha(5); a6=ha(6);
%%
alphaxx = 0:0.1:0.8;

%%  mult
alpha2 = alphaxx;  % fraction of rate of susbtrate transferred to neigbhour grid cells 
mult_path='Mult_transport\output\';
load([mult_path,'Ph_mult_noCorr_',num2str(alpha2(1)),'.mat'])
p1= plot(ax1,Ch(:,1).*1e-12./(rho*Vmic_cc), (1-Y).*dec0./Ch(:,2),'-',...
    'DisplayName','Homogeneous','linewidth',lw, 'Color', clr);hold(ax1, 'on')


warning('off', 'all')
for  id=1:length(alpha2)
load([mult_path,'Ph_mult_noCorr_',num2str(alpha2(id)),'.mat'])
plot(ax1,Cs_mean.*1e-12./(rho*Vmic_cc), (1-Y).*dec_m./Cb_mean,'--',...
    'DisplayName',sprintf('alpha=%1.2f',alpha2(id)),'linewidth',lw, 'Color', LC(id,:));
end
warning('on', 'all')


load([mult_path,'Ph_mult_fullHetero_noCorr_',num2str(alpha2(1)),'.mat'])
plot(a2,Ch(:,1).*1e-12./(rho*Vmic_cc), (1-Y).*dec0./Ch(:,2),'-','linewidth',lw, 'Color', clr);hold(a2, 'on')


warning('off', 'all')
for  id=1:length(alpha2)
load([mult_path,'Ph_mult_fullHetero_noCorr_',num2str(alpha2(id)),'.mat'])
plot(a2,Cs_mean.*1e-12./(rho*Vmic_cc), (1-Y).*dec_m./Cb_mean,'--','linewidth',lw, 'Color', LC(id,:));
end
warning('on', 'all')

%%  MM
alpha2 = alphaxx;  % fraction of rate of susbtrate transferred to neigbhour grid cells 
mm_path='MM_transient\results\output_transport\varing_alpha\';
load([mm_path,'results_v3_Ph_mm_noCorr_',num2str(alpha2(1)),'.mat'])
plot(a3,Ch(:,1).*1e-12./(rho*Vmic_cc), (1-Y).*dec0./Ch(:,2),'-','linewidth',lw, 'Color', clr);hold(a3, 'on')


warning('off', 'all')
for  id=1:length(alpha2)
load([mm_path,'results_v3_Ph_mm_noCorr_',num2str(alpha2(id)),'.mat'])
plot(a3,Cs_mean.*1e-12./(rho*Vmic_cc), (1-Y).*dec_m./Cb_mean,'--','linewidth',lw, 'Color', LC(id,:));
end
warning('on', 'all')


load([mm_path,'results_v3_Ph_mm_fullHetero_noCorr_',num2str(alpha2(1)),'.mat'])
plot(a4,Ch(:,1).*1e-12./(rho*Vmic_cc), (1-Y).*dec0./Ch(:,2),'-','linewidth',lw, 'Color', clr);hold(a4, 'on')


warning('off', 'all')
for  id=1:length(alpha2)
load([mm_path,'results_v3_Ph_mm_fullHetero_noCorr_',num2str(alpha2(id)),'.mat'])
plot(a4,Cs_mean.*1e-12./(rho*Vmic_cc), (1-Y).*dec_m./Cb_mean,'--','linewidth',lw, 'Color', LC(id,:));
end
warning('on', 'all')


%% Inv MM

alpha2 = alphaxx;  % fraction of rate of susbtrate transferred to neigbhour grid cells 
inv_path='MM_Inv_transient\output\';
load([inv_path,'Ph_mm_inv_noCorr_',num2str(alpha2(1)),'.mat'])
plot(a5,Ch2(:,1).*1e-12./(rho*Vmic_cc), (1-Y).*decmminv./Ch2(:,2),'-','linewidth',lw, 'Color', clr);hold(a5, 'on')


warning('off', 'all')
for  id=1:length(alpha2)
load([inv_path,'Ph_mm_inv_noCorr_',num2str(alpha2(id)),'.mat'])
plot(a5,Cs_mean.*1e-12./(rho*Vmic_cc), (1-Y).*dec_m./Cb_mean,'--','linewidth',lw, 'Color', LC(id,:));
end
warning('on', 'all')


load([inv_path,'Ph_mm_inv_fullhetero_noCorr_',num2str(alpha2(1)),'.mat'])
plot(a6,Ch2(:,1).*1e-12./(rho*Vmic_cc), (1-Y).*decmminv./Ch2(:,2),'-','linewidth',lw, 'Color', clr);hold(a6, 'on')


warning('off', 'all')
for  id=1:length(alpha2)
load([inv_path,'Ph_mm_inv_fullhetero_noCorr_',num2str(alpha2(id)),'.mat'])
plot(a6,Cs_mean.*1e-12./(rho*Vmic_cc), (1-Y).*dec_m./Cb_mean,'--',...
    'DisplayName',sprintf('alpha=%1.2f',alpha2(id)),'linewidth',lw, 'Color', LC(id,:));
end
warning('on', 'all')
%%
%%
hxlable=xlabel(a6,'$\overline{C}_S$ (mgC/gSoil)','Interpreter','latex','FontSize', fsz, 'FontName', ftype2);
hx=hxlable.Position;
hxlable.Position= [hx(1)-90,hx(2)-0.006,hx(3)];

hylabel1=ylabel(a3,'$SGR$ \ $\mathrm{(h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2);
hlegend=legend(a2,p1,'Homogeneous');
set(hlegend,'box', 'off','Location', 'northwest', 'Interpreter','latex','FontSize', fsz, 'FontName', ftype1);

% hlegend=legend(ax1,'show');
% set(hlegend,'box', 'off','Location', 'northwest', 'Interpreter','latex','FontSize', fsz, 'FontName', ftype1);
% hlegend.NumColumns=4;
% hlegend.Position = [0.3151    0.0790     0.4328    0.0666];
% 
ylim(ha,[0, 0.04]);
set(ha(1:4),'XTickLabel',[]);
set(ha(2:2:6),'YTickLabel',[]);
set(ha,'FontSize', fsz,'FontName',ftype1, 'LineWidth',  alw,...
    'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on',...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);

% Biophysical heterogeneity \newline
% Full heterogeneity
text(ax1,20,0.035, 'Mult','FontSize', fsz, 'FontName', ftype2 )
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
axes(ax1)
c2=colorbar;
set(c2,'Location','southoutside', 'FontName', ftype2, 'FontSize', fsz,...
    'Box','off','TickDirection','out'); 
c2.Position=[ 0.3437    0.1144    0.3744    0.0102];
c2.Limits= [0, max(alpha2)];
set(get(c2,'Title'),'String','alpha','FontName', ftype2, 'FontSize', fsz)
set(get(c2,'Title'),'Position',[70.2000 -30 0])

align_Ylabels(gcf)

set(gcf, 'Color', 'w')
export_fig('Fig9.pdf', '-r300');

