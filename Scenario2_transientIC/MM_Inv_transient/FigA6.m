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

phi = (1+sqrt(5))/2; %godlen ratio
width = 4;     % Width in inches
height = width*phi;    % Height in inches
alw = 0.5;    % AxesLineWidth
fsz = 11;      % FontSize
fszax=fsz-1;
fsz_title= fsz+1;
fsz_ax_lgnd = fsz;
lws=2;
lw = 0.75;      % LineWidth
lwdot= 0.5;
lwdotmfa= 1.5;

msz = 4;       % MarkerSize
clr = [0.1,0.1,0.1];     % Line color
ftype1= 'Helvetica';  %hLegend, gca
ftype2='Times New Roman'; % hTitle, hXLabel, hYLabel
xycolor= [0.1 0.1 0.1]*2;
fig=figure;
pos = get(gcf, 'Position');
set(gcf, 'Position', [20 50 width*100, height*100]); %<- Set size
% [ha, pos] = tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
[ha, pos] = tight_subplot(3,1, [.075 .05],[.25 .03],[.2 .05]);
ax1=ha(1);a2=ha(2);a3=ha(3);

%%
load('Ph_mm_inv_pos_0.mat');
ksmm_inv=mean2(ksmm_inv);
Km_inv=mean2(Km_inv);
kb=mean2(kb);
vart=cov1;
covt=cov4;
dect=ksmm_inv.*Cs_mean.*Cb_mean./(Cb_mean+Km_inv);
II_terms= vart+covt;
sot= II_terms + dect; 
tmp=1:3:200;
plot(ax1,t./24, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws ); hold(ax1,'on')
plot(ax1,t./24, (1-Y).*dect.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa ); 
plot(ax1,t./24, (1-Y).*II_terms.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws, 'Color', clr);   
plot(ax1,t./24, (1-Y).*decmminv.*1e-12./(rho*Vmic_cc),'-','linewidth',lws ); 
plot(ax1,t(tmp(1:end))./24, (1-Y).*vart(tmp(1:end)).*1e-12./(rho*Vmic_cc),'-o','linewidth',lwdot,...
     'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(ax1,t(tmp(1:end))./24, (1-Y).*covt(tmp(1:end)).*1e-12./(rho*Vmic_cc),'-d','linewidth',lwdot ,...
    'MarkerSize',msz,'MarkerEdgeColor',xycolor); 



load('Ph_mm_inv_neg_0.mat');
ksmm_inv=mean2(ksmm_inv);
Km_inv=mean2(Km_inv);
kb=mean2(kb);
vart=cov1;
covt=cov4;
dect=ksmm_inv.*Cs_mean.*Cb_mean./(Cb_mean+Km_inv);
II_terms= vart+covt;
sot= II_terms + dect; 
tmp=1:3:200;
plot(a2,t./24, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws ); hold(a2,'on')
plot(a2,t./24, (1-Y).*dect.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa ); 
plot(a2,t./24, (1-Y).*II_terms.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws, 'Color', clr);   
plot(a2,t./24, (1-Y).*decmminv.*1e-12./(rho*Vmic_cc),'-','linewidth',lws ); 
plot(a2,t(tmp(1:end))./24, (1-Y).*vart(tmp(1:end)).*1e-12./(rho*Vmic_cc),'-o','linewidth',lwdot,...
     'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a2,t(tmp(1:end))./24, (1-Y).*covt(tmp(1:end)).*1e-12./(rho*Vmic_cc),'-d','linewidth',lwdot ,...
    'MarkerSize',msz,'MarkerEdgeColor',xycolor); 

load('Ph_mm_inv_noCorr_0.mat');
ksmm_inv=mean2(ksmm_inv);
Km_inv=mean2(Km_inv);
kb=mean2(kb);
vart=cov1;
covt=cov4;
dect=ksmm_inv.*Cs_mean.*Cb_mean./(Cb_mean+Km_inv);
II_terms= vart+covt;
sot= II_terms + dect; 
tmp=1:3:200;
plot(a3,t./24, (1-Y).*dec_m.*1e-12./(rho*Vmic_cc),'--','linewidth',lws ); hold(a3,'on')
plot(a3,t./24, (1-Y).*dect.*1e-12./(rho*Vmic_cc),':','linewidth',lwdotmfa ); 
plot(a3,t./24, (1-Y).*II_terms.*1e-12./(rho*Vmic_cc),'-.','linewidth',lws, 'Color', clr);   
plot(a3,t./24, (1-Y).*decmminv.*1e-12./(rho*Vmic_cc),'-','linewidth',lws ); 
plot(a3,t(tmp(1:end))./24, (1-Y).*vart(tmp(1:end)).*1e-12./(rho*Vmic_cc),'-o','linewidth',lwdot,...
     'MarkerSize',msz,'MarkerEdgeColor',xycolor); 
plot(a3,t(tmp(1:end))./24, (1-Y).*covt(tmp(1:end)).*1e-12./(rho*Vmic_cc),'-d','linewidth',lwdot ,...
    'MarkerSize',msz,'MarkerEdgeColor',xycolor); 

TT=100;xtick=0:25:TT;
set(ha(1),'XLim',[0,TT],'XTickLabel',[],'YLim',[-0.1,0.15],'YTick',-0.4:0.2:0.2,...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
% axis(ha(1), 'square');
set(ha(2),'XLim',[0,TT],'XTickLabel',[],'YLim',[-0.1,0.15],'YTick',-0.4:0.2:0.2,...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
% axis(ha(2), 'square');
set(ha(3),'XLim',[0,TT],'XTick',xtick,'YLim',[-0.1,0.15],'YTick',-0.4:0.2:0.2,...
    'FontSize', fszax,'FontName',ftype1,...
    'LineWidth',  alw,'Box', 'off',...
    'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
% axis(ha(3), 'square');

% xlabel(a1,'Time (day)','FontSize', fsz, 'FontName', ftype2,'Interpreter','latex')
% xlabel(a2,'Time (day)','FontSize', fsz, 'FontName', ftype2,'Interpreter','latex')
xlabel(a3,'Time (day)','FontSize', fsz, 'FontName', ftype2,'Interpreter','latex')
ylabel(ax1,'$\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a2,'$\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
ylabel(a3,'$\mathrm{(mgC/gSoil \ h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)

align_Ylabels(gcf)

title(ax1,'Positive','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
title(a2,'Negative','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
title(a3,'Uncorrelated','Interpreter','latex','FontSize', fsz, 'FontName', ftype2)
axes(a3);
[hlegend,hobj]=columnlegend(2,{'$\overline{R}_{het}$', 'MFA','$\sum$ SOT',...
    '$\overline{R}_{hom}$','Var term','Cov term'}, 'location','northeast','Interpreter','latex');
set(hlegend,'box', 'off','Location', 'northeast', 'Interpreter','latex','FontSize', fsz_ax_lgnd, 'FontName', ftype1);
objhl = findobj(hobj, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 6); %// set marker size as desired
lpos=hlegend.Position;
lpos(3)=lpos(3)+0.25;
lpos=hlegend.Position;

pos1=get(a3,'Position');
legend_x=( pos1(1)+pos1(3) -pos1(1) -lpos(3))/2  +pos1(1) ;
legend_h= pos1(4)-pos1(2)/1.25;
set(hlegend, 'Position',[legend_x-0.12 legend_h   lpos(3)+0.25    lpos(4)])

str={'a','b','c'};
for i=1:3
    axes(ha(i)); %set the current axes to axes2
    text(TT-10, ha(i).YLim(2) - ha(i).YLim(2)/5, ['(',str{i},')'], ...
        'FontWeight','normal','FontSize', fsz, 'FontName', ftype2)
end
set(gcf, 'Color','w');

%save
export_fig(gcf, 'FigA6.pdf', '-r300');

%%
rmpath(selpath)
rmpath(genpath([selpath,'\Third_party_scripts\']))
rmpath([selpath,'\Scenario2_transientIC\MM_Inv_transient\output\'])
addpath([selpath,'\Scenario2_transientIC\'])


