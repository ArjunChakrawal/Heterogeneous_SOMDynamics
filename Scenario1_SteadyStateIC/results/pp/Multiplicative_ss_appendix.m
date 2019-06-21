
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
rho_wood=1.1;%g/cm3  max carbon density in form of wood 
nx=100;ny=100; %One pixel = 25micron
Vol_mic = 50*50*50; %um^3
Vmic_cc=Vol_mic*1e-12; %cm3
%% homo mult
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
width =7;     % Width in inches
height = 8;    % Height in inches
alw = 0.5;    % AxesLineWidth
fsz = 11;      % FontSize
fszax=fsz-1;
fsz_title= fsz+1;
fsz_ax_lgnd = fsz;
lws=2;
lw = 0.75;      % LineWidth
lwdot= 0.5;
lwdotmfa= 1.5;
msz = 3;       % MarkerSize
clr = [0.1,0.1,0.1];     % Line color
ftype1= 'Helvetica';  %hLegend, gca
ftype2='Times New Roman'; % hTitle, hXLabel, hYLabel
xycolor= [0.1 0.1 0.1]*2;
fig=figure;
pos = get(gcf, 'Position');
set(gcf, 'Position', [20 100 width*100, height*100]); %<- Set size
% [ha, pos] = tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
[ha, pos] = tight_subplot(3,2, [.05 .05],[.225 .05],[.2 .05]);
a1=ha(1);a2=ha(2);a3=ha(3);
a4=ha(4);a5=ha(5);a6=ha(6);

%% Only physical heterogeneity 
load('Ph_mult_pos_ss_1.mat')
plot(a1,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws); hold(a1,'on'); 
plot(a1,t./24./365, (1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa); 
total_cov = (1-Y).*mean2(ks).*cov_cscb.*1e-12./(rho*Vmic_cc); 
temp= total_cov + (1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc);
plot(a1,t./24./365, total_cov,'-.','linewidth',lws);  
plot(a1,t./24./365, (1-Y).*dec0.*1e-12./(rho*Vmic_cc),'-', 'linewidth',lws);  
xlim(a1,[0, 15]);

load('Ph_mult_neg_ss_1.mat')
plot(a3,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws); hold(a3,'on'); 
plot(a3,t./24./365, (1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa); 
total_cov = (1-Y).*mean2(ks).*cov_cscb.*1e-12./(rho*Vmic_cc); 
temp= total_cov + (1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc);
plot(a3,t./24./365, total_cov,'-.','linewidth',lws);  
plot(a3,t./24./365, (1-Y).*dec0.*1e-12./(rho*Vmic_cc),'-', 'linewidth',lws);  
xlim(a3,[0, 15]);

load('Ph_mult_zero_ss_1.mat')
plot(a5,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws); hold(a5,'on'); 
plot(a5,t./24./365, (1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa); 
total_cov = (1-Y).*mean2(ks).*cov_cscb.*1e-12./(rho*Vmic_cc); 
temp= total_cov + (1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc);
plot(a5,t./24./365, total_cov,'-.','linewidth',lws);  
plot(a5,t./24./365, (1-Y).*dec0.*1e-12./(rho*Vmic_cc),'-', 'linewidth',lws);  
xlim(a5,[0, 15]);

%%
load('Ph1ch2_mult_pos_ss.mat')
tmp=1:5:length(t);
cov1=(1-Y).*(cov_csks).*Cb_mean;
cov2=(1-Y).*(cov_cbks).*Cs_mean;
cov3=(1-Y).*(cov_cscb);
cov4=(1-Y).*(E_cs_cb_ks);

plot(a2,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws); hold(a2,'on'); 
plot(a2,t./24./365, (1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa); 
total_cov = (1-Y).*( cov_csks.*Cb_mean + cov_cbks.*Cs_mean...
                      + cov_cscb.*mean2(ks)+ E_cs_cb_ks); 
temp= total_cov + (1-Y).*mean2(ks).*Cs_mean.*Cb_mean;
plot(a2,t./24./365, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws);  
plot(a2,t./24./365, (1-Y).*dec0.*1e-12./(rho*Vmic_cc),'-', 'linewidth',lws);  
plot(a2,t(tmp)./24./365,cov1(tmp).*1e-12./(rho*Vmic_cc),'-x','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a2,t(tmp)./24./365,cov2(tmp).*1e-12./(rho*Vmic_cc),'-o','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a2,t(tmp)./24./365,cov3(tmp).*mean2(ks).*1e-12./(rho*Vmic_cc),'-s','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a2,t(tmp)./24./365,cov4(tmp).*1e-12./(rho*Vmic_cc),'-d','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
xlim(a2,[0, 15]);


load('Ph1ch2_mult_neg_ss.mat')
cov1=(1-Y).*(cov_csks).*Cb_mean;
cov2=(1-Y).*(cov_cbks).*Cs_mean;
cov3=(1-Y).*(cov_cscb);
cov4=(1-Y).*(E_cs_cb_ks);

plot(a4,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws); hold(a4,'on'); 
plot(a4,t./24./365, (1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa); 
total_cov = (1-Y).*( cov_csks.*Cb_mean + cov_cbks.*Cs_mean...
                      + cov_cscb.*mean2(ks)+ E_cs_cb_ks); 
temp= total_cov + (1-Y).*mean2(ks).*Cs_mean.*Cb_mean;
plot(a4,t./24./365, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws);  
plot(a4,t./24./365, (1-Y).*dec0.*1e-12./(rho*Vmic_cc),'-', 'linewidth',lws);  
plot(a4,t(tmp)./24./365,cov1(tmp).*1e-12./(rho*Vmic_cc),'-x','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a4,t(tmp)./24./365,cov2(tmp).*1e-12./(rho*Vmic_cc),'-o','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a4,t(tmp)./24./365,cov3(tmp).*mean2(ks).*1e-12./(rho*Vmic_cc),'-s','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a4,t(tmp)./24./365,cov4(tmp).*1e-12./(rho*Vmic_cc),'-d','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
xlim(a4,[0, 15]);

load('Ph1ch2_mult_zero_ss.mat')
cov1=(1-Y).*(cov_csks).*Cb_mean;
cov2=(1-Y).*(cov_cbks).*Cs_mean;
cov3=(1-Y).*(cov_cscb);
cov4=(1-Y).*(E_cs_cb_ks);

plot(a6,t./24./365, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws); hold(a6,'on'); 
plot(a6,t./24./365, (1-Y).*mean2(ks).*Cs_mean.*Cb_mean.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa); 
total_cov = (1-Y).*( cov_csks.*Cb_mean + cov_cbks.*Cs_mean...
                      + cov_cscb.*mean2(ks)+ E_cs_cb_ks); 
temp= total_cov + (1-Y).*mean2(ks).*Cs_mean.*Cb_mean;
plot(a6,t./24./365, total_cov.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws);  
plot(a6,t./24./365, (1-Y).*dec0.*1e-12./(rho*Vmic_cc),'-', 'linewidth',lws);  
plot(a6,t(tmp)./24./365,cov1(tmp).*1e-12./(rho*Vmic_cc),'-x','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a6,t(tmp)./24./365,cov2(tmp).*1e-12./(rho*Vmic_cc),'-o','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a6,t(tmp)./24./365,cov3(tmp).*mean2(ks).*1e-12./(rho*Vmic_cc),'-s','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a6,t(tmp)./24./365,cov4(tmp).*1e-12./(rho*Vmic_cc),'-d','linewidth',lwdot,'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
xlim(a6,[0, 15]);

TT=10;
set(ha(1),'XLim',[0,TT],'XTickLabel',[],...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
% axis(ha(1), 'square');
set(ha(2),'XLim',[0,TT],'XTickLabel',[],...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
% axis(ha(2), 'square');
set(ha(3),'XLim',[0,TT],'XTickLabel',[],...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);

set(ha(4),'XLim',[0,TT],'XTickLabel',[],...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
% axis(ha(1), 'square');
set(ha(5),'XLim',[0,TT],...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
% axis(ha(2), 'square');
set(ha(6),'XLim',[0,TT],...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
% axis(ha(3), 'square');


hxlable=xlabel(a5,'Time (year)','FontSize', fsz, 'FontName', ftype2,'Interpreter','latex');
hx=hxlable.Position;
hxlable.Position=[hx(1)+5.5,hx(2)-0.00001,hx(3)];

ylabel(a1,'$\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a3,'$\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a5,'$\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
align_Ylabels(gcf)

axes(a1);
t1=text(-2.5,1e-4,'Positive','Interpreter','latex','FontSize', fsz, 'FontName', ftype2, 'FontWeight',   'bold' );
t1.Rotation=90;

axes(a3);
t1=text(-2.5,1e-4,'Negative','Interpreter','latex','FontSize', fsz, 'FontName', ftype2, 'FontWeight',   'bold' );
t1.Rotation=90;
axes(a5);
t1=text(-2.5,1e-4,'Uncorrelated','Interpreter','latex','FontSize', fsz, 'FontName', ftype2, 'FontWeight',   'bold' );
t1.Rotation=90;

ta1=title(a1, 'Biophysical heterogeneity','FontSize', fsz_title, 'FontName', ftype2,'FontWeight' , 'normal' );
tp1=get(ta1, 'Position');
set(ta1, 'Position', [tp1(1), 8.75e-4])
ta2=title(a2, 'Full heterogeneity','FontSize', fsz_title, 'FontName', ftype2,'FontWeight' , 'normal');
tp2=get(ta2, 'Position');
set(ta2, 'Position', [tp2(1),8.75e-4])


axes(a6); %set the current axes to axes2
[cl,hobj]=columnlegend(2,{'$\overline{R}_{het}$', 'MFA','$ \sum $ HOT','$\overline{R}_{hom}$',...
    '$(1-Y)\overline{C}_b$ $\overline{k_{s,mult}^{\prime} C_s^{\prime}}$',...
    '$(1-Y)\overline{C}_s$ $\overline{k_{s,mult}^{\prime} C_b^{\prime}}$',...
    '$(1-Y)\overline{k}_{s,mult}$ $\overline{C_s^{\prime} C_b^{\prime}}$',...
    '$(1-Y)\overline{k_{s,mult}^{\prime} C_s^{\prime} C_b^{\prime}}$'}, 'location','northeast','Interpreter','latex');
objhl = findobj(hobj, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 6); %// set marker size as desired
set(cl,'box', 'off','FontSize', fsz_ax_lgnd,'FontName', ftype2 );
lpos=cl.Position;
set(cl, 'Position',[lpos(1)-0.15 lpos(2)-0.275   lpos(3)    lpos(4)])

set(gcf, 'Color','w');
str={'a','b','c'};k=1;
for i=1:2:6
    axes(ha(i)); %set the current axes to axes2
    text(TT-1, ha(i).YLim(2) - ha(i).YLim(2)/2.5, ['(',str{k},')'], ...
        'FontWeight','normal','FontSize', fsz, 'FontName', ftype2)
    k=k+1;
end
str={'d','e','f'};
k=1;
for i=2:2:6
    
    axes(ha(i)); %set the current axes to axes2
    text(TT-1, ha(i).YLim(2) - ha(i).YLim(2)/2.5, ['(',str{k},')'], ...
        'FontWeight','normal','FontSize', fsz, 'FontName', ftype2)
    k=k+1;
end

export_fig(gcf,'Mult_SOT_supplimentary.pdf','-r300');
%%
rmpath(selpath)
rmpath([selpath,'\Scenario1_SteadyStateIC\results\'])