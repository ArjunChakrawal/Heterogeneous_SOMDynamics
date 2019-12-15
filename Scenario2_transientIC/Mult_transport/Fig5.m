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

phi = (1 + sqrt(5)) / 2; %godlen ratio
width = 4; % Width in inches
height = width * phi; % Height in inches
alw = 0.5; % AxesLineWidth
fsz = 10; % FontSize
fsz_title = fsz + 1;
fsz_ax_lgnd = fsz;
lw = 1.0; % LineWidth
msz = 8; % MarkerSize
clr = 'k'; % Line color
ftype1 = 'Times New Roman'; %hLegend, gca
ftype2 = 'Times New Roman'; % hTitle, hXLabel, hYLabel
xycolor = [0.2, 0.2, 0.2];
fig = figure;
pos = get(gcf, 'Position');
set(gcf, 'Position', [100, 20, width * 100, height * 100]); %<- Set size
% [ha, pos] = tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
[ha, pos] = tight_subplot(4, 2, [.035, .05], [.15, .1], [.2, .05]);

ax1 = ha(1);
a2 = ha(3);
a3 = ha(5);
a4 = ha(7);
a6 = ha(2);
a7 = ha(4);
a8 = ha(6);
a9 = ha(8);
% a1=ha(1);a2=ha(3);a3=ha(5);
% a6=ha(2);a7=ha(4);a8=ha(6);
TT = 100;
xtick = 0:25:TT;

%% Only physical heterogeneity
load('output\Ph_mult_pos_0.mat')
km = mean2(km);
ksmm = mean2(ksmm);
plot(ax1, ttime./24, Cs_mean.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(ax1, 'on');
plot(a2, ttime./24, Cb_mean.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a2, 'on');
plot(a3, ttime./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a3, 'on');
plot(a4,ttime./24, total_cov.*1e-12./(rho * Vmic_cc),'--','linewidth',lw, 'Color', clr); hold(a4,'on');   


load('output\Ph_mult_neg_0.mat')
km = mean2(km);
ksmm = mean2(ksmm);
plot(ax1, ttime./24, Cs_mean.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a2, ttime./24, Cb_mean.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a3, ttime./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a4,ttime./24, total_cov.*1e-12./(rho * Vmic_cc),'-.','linewidth',lw, 'Color', clr); hold(a4,'on');   

load('output\Ph_mult_noCorr_0.mat')
km = mean2(km);
ksmm = mean2(ksmm);
plot(ax1, ttime./24, Cs_mean.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a2, ttime./24, Cb_mean.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a3, ttime./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a4,ttime./24, total_cov.*1e-12./(rho * Vmic_cc),':','linewidth',lw, 'Color', clr); hold(a4,'on');   

% homo
plot(ax1, ttime./24, Chh(:, 1).*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
plot(a2, ttime./24, Chh(:, 2).*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
plot(a3, ttime./24, (1 - Y).*dec0h.*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
% plot(a5,ttime./24, Ch(:,3).*1e-12./(rho*Vmic_cc),'-','linewidth',lw, 'Color', clr);

%% physical+ chem heterogeneity
load('output\Ph_mult_fullHetero_pos_0.mat')

km = mean2(km);
ksmm = mean2(ksmm);
plot(a6, ttime./24, Cs_mean.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a6, 'on');
plot(a7, ttime./24, Cb_mean.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a7, 'on');
plot(a8, ttime./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a8, 'on');
plot(a9, ttime./24, total_cov.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a9, 'on');


load('output\Ph_mult_fullHetero_neg_0.mat')
km = mean2(km);
ksmm = mean2(ksmm);
plot(a6, ttime./24, Cs_mean.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a7, ttime./24, Cb_mean.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a8, ttime./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a9, ttime./24, total_cov.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);


load('output\Ph_mult_fullHetero_noCorr_0.mat')
km = mean2(km);
ksmm = mean2(ksmm);
plot(a6, ttime./24, Cs_mean.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a7, ttime./24, Cb_mean.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a8, ttime./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a9, ttime./24, total_cov.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);

% homo
plot(a6, ttime./24, Chh(:, 1).*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
plot(a7, ttime./24, Chh(:, 2).*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
plot(a8, ttime./24, (1 - Y).*dec0h.*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
% plot(a10,ttime./24, Ch(:,3).*1e-12./(rho*Vmic_cc),'-','linewidth',lw, 'Color', clr);

%%
xlim(ha, [0, TT]);
ylabel(ax1, '$\overline{C}_S$ (mgC/gSoil)', 'Interpreter', 'latex', 'FontSize', fsz, 'FontName', ftype2)
ylabel(a2, '$\overline{C}_B$ (mgC/gSoil)', 'Interpreter', 'latex', 'FontSize', fsz, 'FontName', ftype2)
ylabel(a3, '$\overline{R}$  $\mathrm{(mgC/gSoil \ h^{-1})}$', 'Interpreter', 'latex', 'FontSize', fsz, 'FontName', ftype2)
ylabel(a4, '$\sum$ SOT', 'Interpreter', 'latex')
% ylabel(a5,'$\overline{CO}_2$ (mgC/gSoil)','Interpreter','latex')

hxlable = xlabel(a4, 'Time (day)', 'FontSize', fsz, 'FontName', ftype2, 'Interpreter', 'latex');
hxlable.Position = [100, -0.4232, -1.0000];

ta1 = title(ax1, 'Biophysical heterogeneity', 'FontSize', fsz_title, 'FontName', ftype2, 'FontWeight', 'normal');
ta2 = title(a6, 'Full heterogeneity', 'FontSize', fsz_title, 'FontName', ftype2, 'FontWeight', 'normal');
ta1.Position = [50.0001, 166.1538, 0];
ta2.Position = ta1.Position;


ylim(ax1, [0, 150]); ylim(a6, [0, 150])
ylim(a2, [0, 45]); ylim(a7, [0, 45])
ylim(a3, [0, 0.2]); ylim(a8, [0, 0.2])
ylim([a4,a9], [-0.15, 0.1]);


set(ax1, 'XTickLabel', [], 'YTick', 0:50:150, 'FontSize', fsz, 'FontName', ftype1, ...
    'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);
set(a2, 'XTickLabel', [], 'YTick', 0:15:50, 'FontSize', fsz, 'FontName', ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);
set(a3, 'XTickLabel', [], 'YTick', 0:0.1:0.2, 'FontSize', fsz, 'FontName', ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);
set(a4, 'XTick', xtick, 'YTick', -0.15:0.1:0.15, 'FontSize', fsz, 'FontName', ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);
set(a6, 'XTickLabel', [], 'YTickLabel', [], 'YTick', 0:50:150, 'FontSize', fsz, 'FontName', ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);
set(a7, 'XTickLabel', [], 'YTickLabel', [], 'YTick', 0:15:45, 'FontSize', fsz, 'FontName', ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);
set(a8, 'XTickLabel', [], 'YTickLabel', [], 'YTick', 0:0.1:0.2, 'FontSize', fsz, 'FontName', ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);
set(a9, 'XTick', xtick, 'YTickLabel', [], 'YTick', -0.15:0.1:0.15, 'FontSize', fsz, 'FontName', ftype1, 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);

set(gcf, 'Color', 'w');
numbers = 1:26;
lettimeers = lower(char(numbers+64));
for i = 1:length(ha)
    axes(ha(i)); %set the current axes to axes2
    text(TT-10, ha(i).YLim(2)-ha(i).YLim(2)/10, ['(', lettimeers(i), ')'], 'FontWeight', 'normal', 'FontSize', fsz, 'FontName', ftype2)
end

align_Ylabels(gcf)


hlegend = legend(ax1, 'Positive', 'Negative', 'Uncorrelated', 'Homogeneous');
hlegend.NumColumns = 2;
set(hlegend, 'box', 'off', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', fsz_ax_lgnd, 'FontName', ftype1);
set(hlegend, 'Position', [0.3130, 0.0149, 0.3874, 0.0406])

export_fig('figure5.pdf', '-r300');
%%
rmpath(selpath)
rmpath(genpath([selpath,'\Third_party_scripts\']))
rmpath(genpath([selpath,'\Scenario2_transientIC\Mult_transport\']))