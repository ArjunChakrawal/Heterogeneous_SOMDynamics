clc
% close all
clearvars;

rho= 1.65; %g/cm3 soil bulk density 
Vol_mic = 50*50*50; %um^3
Vmic_cc=Vol_mic*1e-12; %cm3
selpath = uigetdir;
addpath(selpath)
addpath(genpath([selpath,'\Third_party_scripts\']))
addpath([selpath,'\Scenario1_SteadyStateIC\results\'])

load('Ph_mult_pos_ss_1.mat')
spCb= cb.*1e-12./(rho*Vmic_cc);
csp= cs.*1e-12./(rho*Vmic_cc);
load('Ph_mult_neg_ss_1.mat')
csn= cs.*1e-12./(rho*Vmic_cc);
load('Ph_mult_zero_ss_1.mat')
csz= cs.*1e-12./(rho*Vmic_cc);

homo= ones(100,100).*(mean2(kb)/(Y*mean2(ks))).*1e-12./(rho*Vmic_cc);
%% Figure properties
width = 8;     % Width in inches
height = 3;    % Height in inches
alw = 0.01;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
clr = 'k';     % Line color
fonttype = 'Times New Roman'; 
cxl = 355;
bxl = 0.12;
fig=figure;
fig.Name= 'maps';
pos = get(gcf, 'Position');
set(gcf, 'Position', [20 100 width*100, height*100]); %<- Set size

xtick= [];
ytick= [];
%%
a1=subplot(1,5,1);
surf(a1,homo);  shading interp; view(0,90); 
c1=colorbar; 
axis square; axis off

a2=subplot(1,5,2);
[~,ch] = contourf(csp); ch.LineStyle = 'none'; ch.LevelStep = ch.LevelStep/10;
c2=colorbar; 
axis off
axis square

a3=subplot(1,5,3);
[~,ch] = contourf(csn); ch.LineStyle = 'none'; ch.LevelStep = ch.LevelStep/10;
c3=colorbar; 
axis off
axis square

a4=subplot(1,5,4);
[~,ch] = contourf(csz); ch.LineStyle = 'none'; ch.LevelStep = ch.LevelStep/50;
c4=colorbar; 
axis off
axis square

a5=subplot(1,5,5);
[~,ch] = contourf(spCb); ch.LineStyle = 'none'; ch.LevelStep = ch.LevelStep/50;
c5=colorbar; 
axis off
axis square


%% labels and ticks
t1=title(a1, 'Homogeneous','FontWeight','Normal'); 
t2=title(a2, 'Positive','FontWeight','Normal');
t3=title(a3, 'Negative','FontWeight','Normal');
t4=title(a4, 'Uncorrelated','FontWeight','Normal');
t5=title(a5, 'Microbial C','FontWeight','Normal');


set(a1,'XTick',xtick,'YTick',ytick,'FontSize', fsz,'FontName',fonttype, 'LineWidth', alw);
set(a2,'XTick',xtick,'YTick',ytick,'FontSize', fsz,'FontName',fonttype, 'LineWidth', alw);
set(a3,'XTick',xtick,'YTick',ytick,'FontSize', fsz,'FontName',fonttype, 'LineWidth', alw);
set(a4,'XTick',xtick,'YTick',ytick,'FontSize', fsz,'FontName',fonttype, 'LineWidth', alw);
set(a5,'XTick',xtick,'YTick',ytick,'FontSize', fsz,'FontName',fonttype, 'LineWidth', alw);

set(c1,'Location','southoutside', 'FontName', fonttype, 'FontSize', fsz,...
    'Ticks',mean2(homo),'TickLabels',{sprintf('%1.1f',mean2(homo))},...
    'Box','off','TickDirection','out'); 
set(c2, 'Location','southoutside', 'FontName', fonttype, 'FontSize', fsz,'Box','off', ...
    'Ticks', linspace(min(csp(:)),10.5,3),'TickLabels',{sprintf('%1.1f',min(csp(:))),'6.7','10.6'},'TickDirection','out');

set(c3, 'Location','southoutside', 'FontName', fonttype, 'FontSize', fsz, 'Box','off',...
    'Ticks', linspace(min(csn(:)),7.6,3),'TickLabels',{sprintf('%1.1f',min(csn(:))),'5.6','7.6'},'TickDirection','out');
set(c4, 'Location','southoutside', 'FontName', fonttype, 'FontSize', fsz, 'Box','off',...
    'Ticks', linspace(min(csz(:)),23,3),'TickLabels',{sprintf('%1.0f',min(csz(:))),'11','23.5'},'TickDirection','out');
set(c5, 'Location','southoutside', 'FontName', fonttype, 'FontSize', fsz, 'Box','off',...
    'Ticks', linspace(min(spCb(:)),2.4,3),'TickLabels',{'1E-4','1.2','2.4'},...
    'TickDirection','out');

colormap jet

colormap(a1,[0 1 0.9])

colormap jet

text(a1,40, 130,'(a)', 'FontName', fonttype, 'FontSize', fsz)
text(a2,40, 130,'(b)', 'FontName', fonttype, 'FontSize', fsz)
text(a3,40, 130,'(c)', 'FontName', fonttype, 'FontSize', fsz)
text(a4,40, 130,'(d)', 'FontName', fonttype, 'FontSize', fsz)
text(a5,40, 130,'(e)', 'FontName', fonttype, 'FontSize', fsz)
text(a3,5, -60,'(mgC/gSoil)','Interpreter','latex','FontName', fonttype, 'FontSize', fsz)

set(gcf,'Color','w');
set(gcf, 'Renderer', 'opengl');
export_fig(gcf,'maps_ss.pdf','-r300');

%%
rmpath(selpath)
rmpath([selpath,'\Scenario1_SteadyStateIC\results\'])
rmpath(genpath([selpath,'\Third_party_scripts\']))