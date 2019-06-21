
%   in inches% Figure properties
clear all
close all
selpath = uigetdir;
addpath(selpath)
addpath([selpath,'\Scenario1_SteadyStateIC\Spatial_field'])
addpath([selpath,'\Scenario2_transientIC\Spatial_field'])

phi= (1+sqrt(5))/2;
height = 3.5;     %   in inches
width = 7;
alw = 0.5;    % AxesLineWidth
fsz = 11;      % FontSize
fszax=fsz-1;

fsz_title= fsz+1;
fsz_ax_lgnd = fsz;
lw = 1.0;      % LineWidth
msz = 8;       % MarkerSize
clr = 'k';     % Line color
ftype1= 'Helvetica';  %hLegend, gca
ftype2='Times New Roman'; % hTitle, hXLabel, hYLabel
xycolor= [0.2 0.2 0.2];


fig=figure;
fig.Name = 'Ksmulthist';
pos = get(gcf, 'Position');
set(gcf, 'Position', [20 100 width*100, height*100]); %<- Set size

load('ksm1.mat')
[n, ed]=histcounts(ksm1, 20);
b_center=ed(1:end-1)+diff(ed)/2;
d=diff(ed);
a1=subplot(1,2,1);
H1= area(a1,b_center, n/(sum(n)*d(1)));   hold on
H1.FaceColor= [1,1,1]*0.2;H1.EdgeColor='none';
load('ksm2.mat');
[n, ed]=histcounts(ksm2, 20);
b_center=ed(1:end-1)+diff(ed)/2;
d=diff(ed);
H2=area(a1,b_center, n/(sum(n)*d(1)));   hold on
H2.FaceColor=[1,1,1]*0.6;H2.EdgeColor='none';

plot(ones(1,25).*mean2(ksm1),linspace(0, max(n/(sum(n)*d(1))), 25), '--k', 'LineWidth', lw)

set(a1,'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on','XColor', xycolor, 'YColor', xycolor);
axis square
xlabel(a1,'$k_{s,mult}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a1,'Normalised Rel. freq.','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)

%%
a2=subplot(1,2,2);
load('ksmmh1.mat')
[n, ed]=histcounts(ksmmh1, 20);
b_center=ed(1:end-1)+diff(ed)/2;
d=diff(ed);
H1= area(a2,b_center, n/(sum(n)*d(1)));   hold on
H1.FaceColor= [1,1,1]*0.2;H1.EdgeColor='none';
set(a2,'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on','XColor', xycolor, 'YColor', xycolor);
axis square
plot(a2,ones(1,25).*mean2(ksmmh1),linspace(0, max(n/(sum(n)*d(1))), 25), '--k', 'LineWidth', lw)
xlabel(a2,'$k_{s,mm}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)


set(gcf,'Color','w');

export_fig(gcf,'Ks_hist2.pdf','-r300');



%%
rmpath(selpath)
rmpath([selpath,'\Scenario1_SteadyStateIC\Spatial_field'])
rmpath([selpath,'\Scenario2_transientIC\Spatial_field'])
