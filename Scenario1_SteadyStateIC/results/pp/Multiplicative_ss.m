clc
close all
clearvars;

selpath = uigetdir;
addpath(selpath)
addpath([selpath,'\Scenario1_SteadyStateIC\results\'])
ksmm=0.018;
kb=0.00028;  %h-1
Y=0.31;
km=25;
rho= 1.65; %g/cm3 soil bulk density 
rho_wood=1.1;%g/cm3  max carboff density in form of wood 
nx=100;ny=100; %offe pixel = 25microff

Vol_mic = 50*50*50; %um^3
Vmic_cc=Vol_mic*1e-12; %cm3
Vdomain_cc=Vol_mic*nx*ny*1e-12 ; % cm^3

%% homo
load('Ph_mult_neg_ss_1.mat')
cs0=mean2(cs);cb0=mean2(cb);co20=0;
ks=mean2(ks);
kb=mean2(kb);
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(t)).*Imic;
[~,C1]= DEC_multiplicative(cs0,cb0,co20,t, ks, kb, Y,St, options);
dec0=ks.*C1(:,1).*C1(:,2);
%% Figure properties
phi = (1+sqrt(5))/2; %godlen ratio
width = 4;     % Width in inches
height = width*phi;    % Height in inches
alw = 0.5;    % AxesLineWidth
fsz = 10;      % FontSize
fsz_title= fsz+1;
fsz_ax_lgnd = fsz;
lw = 1.0;      % LineWidth
msz = 8;       % MarkerSize
clr = 'k';     % Line color
ftype1= 'Helvetica';  %hLegend, gca
ftype2='Times New Roman'; % hTitle, hXLabel, hYLabel
xycolor= [0.2 0.2 0.2];
fig=figure;
pos = get(gcf, 'Position');
set(gcf, 'Position', [20 100 width*100, height*100]); %<- Set size
% [ha, pos] = tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
[ha, pos] = tight_subplot(4, 2, [.035 .05],[.15 .1],[.2 .05]);
a1=ha(1);a2=ha(3);a3=ha(5);a4=ha(7);
a6=ha(2);a7=ha(4);a8=ha(6);a9=ha(8);

TT=10;
xtick= 0:5:15;
% cstck=[0,40,80];
%% only physical heterogeneity 
load('Ph_mult_pos_ss_1.mat');
plot(a1,t./24./365, Cs_mean.*1e-12./(rho*Vmic_cc),'--','linewidth',lw, 'Color', clr);   hold(a1,'on');  
plot(a2,t./24./365, Cb_mean.*1e-12./(rho*Vmic_cc),'--','linewidth',lw, 'Color', clr);  hold(a2,'on'); 
plot(a3,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lw, 'Color', clr);  hold(a3,'on'); 
plot(a4,t./24./365, (1-Y).*mean2(ks).*cov_cscb.*1e-12./(rho*Vmic_cc),'--','linewidth',lw, 'Color', clr); hold(a4,'on');   

load('Ph_mult_neg_ss_1.mat');
plot(a1,t./24./365, Cs_mean.*1e-12./(rho*Vmic_cc),'-.','linewidth',lw, 'Color', clr);  
plot(a2,t./24./365, Cb_mean.*1e-12./(rho*Vmic_cc),'-.','linewidth',lw, 'Color', clr);  
plot(a3,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'-.','linewidth',lw, 'Color', clr);  
plot(a4,t./24./365, (1-Y).*mean2(ks).*cov_cscb.*1e-12./(rho*Vmic_cc),'-.','linewidth',lw, 'Color', clr);  

load('Ph_mult_zero_ss_1.mat');
plot(a1,t./24./365, Cs_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lw, 'Color', clr);  
plot(a2,t./24./365, Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lw, 'Color', clr);  
plot(a3,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),':','linewidth',lw, 'Color', clr);  
plot(a4,t./24./365, (1-Y).*mean2(ks).*cov_cscb.*1e-12./(rho*Vmic_cc),':','linewidth',lw, 'Color', clr);  

% homo
x=4;
plot(a1,t./24./365, C1(:,1).*1e-12./(rho*Vmic_cc),'-','linewidth',lw, 'Color', clr);  
plot(a2,t./24./365, C1(:,2).*1e-12./(rho*Vmic_cc),'-', 'linewidth',lw, 'Color', clr);  
plot(a3,t./24./365, (1-Y).*dec0.*1e-12./(rho*Vmic_cc),'-', 'linewidth',lw, 'Color', clr);  


%% physical+ chem heterogeneity 
load('Ph1ch2_mult_pos_ss.mat')
plot(a6,t./24./365, Cs_mean.*1e-12./(rho*Vmic_cc),'--','linewidth',lw, 'Color', clr);  hold(a6,'on'); 
plot(a7,t./24./365, Cb_mean.*1e-12./(rho*Vmic_cc),'--','linewidth',lw, 'Color', clr);  hold(a7,'on'); 
plot(a8,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lw, 'Color', clr);  hold(a8,'on'); 
total_cov = (1-Y).*( cov_csks.*Cb_mean + cov_cbks.*Cs_mean...
                      + cov_cscb.*mean2(ks)+ E_cs_cb_ks); 
temp= total_cov + (1-Y).*mean2(ks).*Cs_mean.*Cb_mean;
plot(a9,t./24./365, total_cov.*1e-12./(rho*Vmic_cc),'--','linewidth',lw, 'Color', clr);  hold(a9,'on'); 

load('Ph1ch2_mult_neg_ss.mat')
plot(a6,t./24./365, Cs_mean.*1e-12./(rho*Vmic_cc),'-.','linewidth',lw, 'Color', clr);  
plot(a7,t./24./365, Cb_mean.*1e-12./(rho*Vmic_cc),'-.','linewidth',lw, 'Color', clr);  
plot(a8,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'-.','linewidth',lw, 'Color', clr);  
total_cov = (1-Y).*( cov_csks.*Cb_mean + cov_cbks.*Cs_mean...
                      + cov_cscb.*mean2(ks)+ E_cs_cb_ks); 
% temp= total_cov + (1-Y).*mean2(ks).*Cs_mean.*Cb_mean;
plot(a9,t./24./365, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lw, 'Color', clr);  

load('Ph1ch2_mult_zero_ss.mat')
plot(a6,t./24./365, Cs_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lw, 'Color', clr);  
plot(a7,t./24./365, Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lw, 'Color', clr);  
plot(a8,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),':','linewidth',lw, 'Color', clr);  
total_cov = (1-Y).*( cov_csks.*Cb_mean + cov_cbks.*Cs_mean...
                      + cov_cscb.*mean2(ks)+ E_cs_cb_ks); 
% temp= total_cov + (1-Y).*mean2(ks).*Cs_mean.*Cb_mean;
plot(a9,t./24./365, total_cov.*1e-12./(rho*Vmic_cc),':','linewidth',lw, 'Color', clr);  


% homo
plot(a6,t./24./365, C1(:,1).*1e-12./(rho*Vmic_cc),'-','linewidth',lw, 'Color', clr);  
plot(a7,t./24./365, C1(:,2).*1e-12./(rho*Vmic_cc),'-', 'linewidth',lw, 'Color', clr);  
plot(a8,t./24./365, (1-Y).*dec0.*1e-12./(rho*Vmic_cc),'-', 'linewidth',lw, 'Color', clr);  


%%
rmpath(selpath)
rmpath([selpath,'\Scenario1_SteadyStateIC\results\'])


%%  figure beautification
xlim(a1,[0, TT]);
xlim(a2,[0, TT]);
xlim(a3,[0, TT]);
xlim(a4,[0, TT]);
xlim(a6,[0, TT]);
xlim(a7,[0, TT]);
xlim(a8,[0, TT]);
xlim(a9,[0, TT]);


cstick= linspace(5, 7.5,4);
cbtick= linspace(.75, 1.1,4);
rtick= linspace(4.5, 7.5,4).*1e-4;
covtick= linspace(-12, 12.5,4).*1e-5;

ylim([a1,a6], [cstick(1),cstick(end)]);

ylim([a2,a7], [cbtick(1),cbtick(end)]);

ylim([a3,a8], [rtick(1),rtick(end)])
 
ylim([a4,a9], [covtick(1),covtick(end)])
 


set(a1,'XTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw,'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a6,'XTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a2,'XTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a7,'XTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a3,'XTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a8,'XTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a4,'XTick',xtick,'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a9,'XTick',xtick,'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);

set(a1,'YTick',cstick,'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a6,'YTick',cstick,'YTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a2,'YTick',cbtick,'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a7,'YTick',cbtick,'YTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a3,'YTick',rtick,'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a8,'YTick',rtick,'YTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a4,'YTick',covtick,'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a9,'YTick',covtick,'YTickLabel',[],'FontSize', fsz_ax_lgnd,'FontName',ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);

ytickformat([a1,a2,a3,a4],'%1.2f')

ylabel(a1,'$\overline{C}_s$ (mgC/gSoil)','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a2,'$\overline{C}_b$ (mgC/gSoil)','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a3,'$\overline{R}$  $\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a4,'$\sum$ HOT','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)

hxlable=xlabel(a4,'Time (year)','FontSize', fsz, 'FontName', ftype2,'Interpreter','latex');
hx=hxlable.Position;
hxlable.Position=[hx(1)+5.5,hx(2)-1e-5,hx(3)];

ta1=title(a1, 'Biophysical heterogeneity','FontSize', fsz_title, 'FontName', ftype2,'FontWeight' , 'normal' );
tp1=get(ta1, 'Position');
set(ta1, 'Position', [tp1(1), 8])
ta2=title(a6, 'Full heterogeneity','FontSize', fsz_title, 'FontName', ftype2,'FontWeight' , 'normal');
tp2=get(ta2, 'Position');
set(ta2, 'Position', [tp2(1),8])

hlegend=legend(a1,'Positive', 'Negative','Uncorrelated','Homogeneous');
set(hlegend,'box', 'off','Location', 'northeast', 'Interpreter','latex','FontSize', fsz_ax_lgnd, 'FontName', ftype1);
hlegend.NumColumns=2;
lpos=hlegend.Position;
pos1=get(a4,'Position'); pos2=get(a9,'Position');
legend_x=( pos2(1)+pos2(3) -pos1(1) -lpos(3))/2  +pos1(1) ;
legend_h= pos1(4)-pos1(2)/1;
set(hlegend, 'Position',[legend_x legend_h   lpos(3)    lpos(4)])
align_Ylabels(gcf)
set(gcf, 'Color', 'w')

str={'a','e','b','f','c','g','d','h'};
for i=1:8
    axes(ha(i)); %set the current axes to axes2
    text(TT-2, ha(i).YLim(2) - ha(i).YLim(2)/50, ['(',str{i},')'], 'FontWeight','normal','FontSize', fsz+2, 'FontName', ftype2)
end

%% save
export_fig(gcf,'Mult_ss.pdf','-r300');

