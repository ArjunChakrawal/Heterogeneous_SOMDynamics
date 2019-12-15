clc
close all
clearvars;
LC = linspecer(4);

selpath = uigetdir;
addpath(selpath)
addpath(genpath([selpath,'\Third_party_scripts\']))
addpath([selpath,'\Scenario2_transientIC\MM_Inv_transient\output\'])
addpath([selpath,'\Scenario2_transientIC\'])

%% define model parameters
ksmm = 0.018; % h-1
kb = 0.00028; % h-1
Y = 0.31; % microbial carbon use efficiency
km = 25; % mgC/gSoil
rho = 1.65; % g/cm3 soil bulk density
rho_wood = 1.1; % g/cm3  max carbon density in form of wood
Vol_mic = 50 * 50 * 50; % volume of a microsite in µm^3 (each grid cell)
Vmic_cc = Vol_mic * 1e-12; %cm3
km = ((km / 1000) .* rho .* 1e-12 .* 1e15) .* Vol_mic; %fgC

% spatial discretization
Lx = 5000;
Ly = 5000; % um
Nx = 100;
Ny = 100;
N = Nx * Ny;
x = linspace(0, Lx, Nx); % mesh points in x dir
y = linspace(0, Lx, Ny); % mesh points in y dir
dx = x(2) - x(1);
dy = y(2) - y(1);

% temporal discretization
T = 24 * 200; % Simulation time 200 days with daily time steps
dt = 24;
Nt = round(T/dt) + 1;
t = 0:dt:T; % Mesh points in time

% transport parameter
Iteration = 'SOR';
% fraction of rate of susbtrate transferred to neigbhour grid cells
alpha = 0:0.1:0.8;

%% Estimation of kinetic parameter
Vdomain_cc = Vol_mic * Nx * Ny * 1e-12; % cm^3
soil_domain = (rho * Vdomain_cc);
faccs = 12.121212; % initial SOM fraction
total_Cs = soil_domain * faccs * 0.01 * 1e15; % total amount of SOM in domain
total_Cb = total_Cs * 0.01; % 1% of total substrate is microbial C
cb0 = (total_Cb / (Nx * Ny)) * 1e-12 ./ (rho * Vmic_cc);
Imic = (total_Cs / (Nx * Ny)) / 200000;
Km_inv = cb0 * 8;
ksmm_inv = ksmm / 4;
Km_inv = ((Km_inv / 1000) .* rho .* 1e-12 .* 1e15) .* Vol_mic; %fgC

load('Spatial_field\field_Csn.mat');
load('Spatial_field\field_Cb.mat');
cs0 = mean2(spCs);
cb0 = mean2(spCb);
co20 = 0;

c(1, :) = [cs0, cb0, co20];
St = ones(1, length(t)) .* Imic; % vector of external C input at each time step
options = odeset('Stats', 'off', 'AbsTol', 1e-9, 'RelTol', 1e-9); % MATLAB ode options
f = @(tt, c)[-ksmm_inv * c(1) * c(2) / (Km_inv + c(2)) + kb * c(2) + interp1(t, St, tt); ...
    Y * ksmm_inv * c(1) * c(2) / (Km_inv + c(2)) - kb * c(2); ...
    (1 - Y) * ksmm_inv * c(1) * c(2) / (Km_inv + c(2))];
[~, Ch] = ode45(f, t, c(1, :), options);
decmminv = ksmm_inv .* Ch(:, 1) .* Ch(:, 2) ./ (Km_inv + Ch(:, 2));

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
set(gcf, 'Position', [20, 100, width * 100, height * 100]); %<- Set size
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

TT = 100;
xtick = 0:25:TT;

%% only physical heterogeneity
load('Ph_mm_inv_pos_0.mat');
plot(ax1, t./24, Cs_mean.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(ax1, 'on');
plot(a2, t./24, Cb_mean.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a2, 'on');
plot(a3, t./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a3, 'on');
plot(a4, t./24, (1 - Y).*total_cov.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a4, 'on');

load('Ph_mm_inv_neg_0.mat');
plot(ax1, t./24, Cs_mean.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a2, t./24, Cb_mean.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a3, t./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a4, t./24, (1 - Y).*total_cov.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);

load('Ph_mm_inv_noCorr_0.mat');
plot(ax1, t./24, Cs_mean.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a2, t./24, Cb_mean.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a3, t./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a4, t./24, (1 - Y).*total_cov.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);

% homo
x = 4;
plot(ax1, t./24, Ch(:, 1).*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
plot(a2, t./24, Ch(:, 2).*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
plot(a3, t./24, (1 - Y).*decmminv.*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);

%% physical+ chem heterogeneity
load('Ph_mm_inv_fullhetero_pos_0.mat')
plot(a6, t./24, Cs_mean.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a6, 'on');
plot(a7, t./24, Cb_mean.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a7, 'on');
plot(a8, t./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a8, 'on');
plot(a9, t./24, total_cov.*1e-12./(rho * Vmic_cc), '--', 'linewidth', lw, 'Color', clr);
hold(a9, 'on');

load('Ph_mm_inv_fullhetero_neg_0.mat')
plot(a6, t./24, Cs_mean.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a7, t./24, Cb_mean.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a8, t./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);
plot(a9, t./24, total_cov.*1e-12./(rho * Vmic_cc), '-.', 'linewidth', lw, 'Color', clr);

load('Ph_mm_inv_fullhetero_noCorr_0.mat')
plot(a6, t./24, Cs_mean.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a7, t./24, Cb_mean.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a8, t./24, (1 - Y).*dec_m.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);
plot(a9, t./24, total_cov.*1e-12./(rho * Vmic_cc), ':', 'linewidth', lw, 'Color', clr);


% homo
plot(a6, t./24, Ch(:, 1).*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
plot(a7, t./24, Ch(:, 2).*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);
plot(a8, t./24, (1 - Y).*decmminv.*1e-12./(rho * Vmic_cc), '-', 'linewidth', lw, 'Color', clr);

%%  figure beautification
for i = 1:8
    xlim(ha(i), [0, TT]);
end
set(ha(3:4), 'Ylim', [0,45],'YTick', 0:15:50);

set(ha(2:2:6), 'YTickLabel', [], 'FontSize', fsz, 'FontName', ftype1, ...
    'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', ...
    'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);

set(ha(1:6), 'XTickLabel', [], 'FontSize', fsz, 'FontName', ftype1, ...
    'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', ...
    'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);

set(ha(7:8), 'XTick', xtick,  'FontSize', fsz, 'FontName', ftype1, ...
    'Ylim', [-0.15, 0.01], ...
    'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', ...
    'TickLength', [.02, .02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XColor', xycolor, 'YColor', xycolor);

set(ha(8), 'YTickLabel', []);
set(ha(5:6), 'Ylim', [0,0.2]);

% ytickformat([ax1,a2,a3,a4],'%1.2f')

ylabel(ax1, '$\overline{C}_S$ (mgC/gSoil)', 'Interpreter', 'latex', 'FontSize', fsz, 'FontName', ftype2)
ylabel(a2, '$\overline{C}_B$ (mgC/gSoil)', 'Interpreter', 'latex', 'FontSize', fsz, 'FontName', ftype2)
ylabel(a3, '$\overline{R}$  $\mathrm{(mgC/gSoil \ h^{-1})}$', 'Interpreter', 'latex', 'FontSize', fsz, 'FontName', ftype2)
ylabel(a4, '$\sum$ SOT', 'Interpreter', 'latex', 'FontSize', fsz, 'FontName', ftype2)

hxlable = xlabel(a4, 'Time (day)', 'FontSize', fsz, 'FontName', ftype2, 'Interpreter', 'latex');
hxlable.Position = [100, -0.1927, -1.0000];
ta1 = title(ax1, 'Biophysical heterogeneity', 'FontSize', fsz_title, 'FontName', ftype2, 'FontWeight', 'normal');
ta2 = title(a6, 'Full heterogeneity', 'FontSize', fsz_title, 'FontName', ftype2, 'FontWeight', 'normal');
ta1.Position=[50.0001  166.1538         0];
ta2.Position=ta1.Position;

hlegend = legend(ax1, 'Positive', 'Negative', 'Uncorrelated', 'Homogeneous');
hlegend.NumColumns = 2;
set(hlegend, 'box', 'off', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', fsz_ax_lgnd, 'FontName', ftype1);
set(hlegend, 'Position', [0.3130, 0.0149, 0.3874, 0.0406])

align_Ylabels(gcf)
set(gcf, 'Color', 'w')

numbers = 1:26;
letters = lower(char(numbers+64));
for i = 1:length(ha)
    axes(ha(i)); %set the current axes to axes2
    text(TT-10, ha(i).YLim(2)-ha(i).YLim(2)/10, ['(', letters(i), ')'], 'FontWeight', 'normal', 'FontSize', fsz, 'FontName', ftype2)
end

%save
export_fig(gcf, 'Fig7.pdf', '-r300');

%%
rmpath(selpath)
rmpath(genpath([selpath,'\Third_party_scripts\']))
rmpath([selpath,'\Scenario2_transientIC\MM_Inv_transient\output\'])
addpath([selpath,'\Scenario2_transientIC\'])
