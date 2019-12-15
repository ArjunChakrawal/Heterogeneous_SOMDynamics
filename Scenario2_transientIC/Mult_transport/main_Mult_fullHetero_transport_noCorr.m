clc;
close all;
clearvars;
selpath = uigetdir;
addpath(selpath)
addpath([selpath,'\Scenario2_transientIC'])

%% define model parameters
ksmm = 0.018; % h-1
kb = 0.00028; % h-1
Y = 0.31; % microbial carbon use efficiency
km = 25; % mgC/gSoil
rho = 1.65; % g/cm3 soil bulk density
rho_wood = 1.1; % g/cm3  max carbon density in form of wood
Vol_mic = 50 * 50 * 50; % volume of a microsite in µm^3 (each grid cell)
Vmic_cc = Vol_mic * 1e-12; %cm3

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


%% %% ksmult estimation and Simulate homogeneous Carbon dynamic for MM kinetics
load('Spatial_field\field_Cs_noCorr.mat');
load('Spatial_field\field_Cb.mat');
cs0=mean2(spCs)*1e-12./(rho*Vmic_cc);cb0=mean2(spCb)*1e-12./(rho*Vmic_cc);co20=0;
Imic=cs0/200000; T=24*1000;
ks=ks_mult_estimation(T,cs0,cb0,co20,kb,ksmm, km, Imic,Y,rho,Vol_mic);

% Homogeneous  case 
cs0=mean2(spCs);cb0=mean2(spCb);co20=0;
Imic=cs0/200000;
options = odeset('Stats','off','AbsTol',1e-6,'RelTol',1e-6);
T=24*200; nstep=200; t=linspace(0,T,nstep);
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(t)).*Imic;

c(1,:)=[cs0,cb0,co20];
f=@(tt,c)[-ks *c(1)*c(2) + kb*c(2) + interp1(t, St,tt);...
    Y*ks *c(1)*c(2) - kb*c(2); (1-Y)*ks*c(1)*c(2)];
[tt, Ch]=ode45(f,t,c(1,:), options);
dec0 = ks .* Ch(:, 1) .* Ch(:, 2) ; %homogeneous rate of decomposition
Cd = Ch .* 1e-12 ./ (rho * Vmic_cc); %Time evolution of carbon mg C at microbsite
figure;
plot(tt./24, Cd(:, 1));
hold on;
plot(tt./24, Cd(:, 2));
plot(tt./24, Cd(:, 3));
xlabel('Time(day)');
ylabel('(mg C/gSoil) ');
legend('Substrate', 'Biomass', 'CO_2'); title('Homogeneous case')
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', ...
    'LineWidth', 0.5, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', [1, 1, 1]*0.3, 'YColor', [1, 1, 1]*0.3);
set(gcf, 'Color', 'w')
figure; plot(t./24, (1 - Y).*dec0.*1e-12./(rho * Vmic_cc))
xlabel('Time(day)');
ylabel('Respirataion rate (mg C/gSoil h^{-1})');
title('Homogeneous case')
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', ...
    'LineWidth', 0.5, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', [1, 1, 1]*0.3, 'YColor', [1, 1, 1]*0.3);
set(gcf, 'Color', 'w')

%% Full heterogeneity only i.e. only substrate and microbial C are spatially heterogeneous
load('Spatial_field\ksm1.mat');

% Relaxation coeff
omega = 0.1;
max_iter = 500;
tol = 1e-6;
abs_tol = 1e-12;
rel_tol = 1e-6;
% Crank–Nicolson theta parameter
theta = 0.5;

% hbar = parfor_progressbar(Nt,'time loop is running, please wait!');  %create the progress bar
mass_balance = zeros(Nt, 1);
tic
for a1 = 1:length(alpha)
    %positive correlation
    ks = ksm1;
    kb = ones(Nx, Ny) .* kb;
    km = ones(Nx, Ny) .* km;
    % The initial condition is easy to fill in:
    us_n = spCs;
    ub_n = spCb;

    SOC = zeros(Nx, Ny); % unknown SOC at new time level
    B = zeros(Nx, Ny);
    co2 = zeros(Nx, Ny);
    uco2_n = zeros(Nx, Ny);
    Cell_SOC = cell(Nt, 1);
    Cell_B = Cell_SOC;
    Cell_co2 = Cell_SOC;

    Cs_mean=zeros(length(t),1);
    Cs_std=Cs_mean;Cs_m3=Cs_mean;
    Cb_std=Cs_mean;
    CO2_mean=Cs_mean;CO2_std=Cs_mean;
    dec_m=Cs_mean;dec_std=Cs_mean;
    cov_cscb=Cs_mean;
    cov_csks=cov_cscb;
    cov_cbks=cov_cscb;
    Cs_m4=zeros(length(t),1);
    Cb_mean=zeros(length(t),1);
    rho_CsCb=zeros(length(t),1);
    rho_Csks= rho_CsCb;
    rho_Cbks=rho_CsCb;
    E_cs_cb_ks=zeros(length(t),1);
    total_cov=zeros(length(t),1);
    %===========INITIAL CONDITION==========================================
    mass_balance(1) = mean2(us_n(:)) + mean2(ub_n(:)) + mean2(uco2_n(:)) ...
        -mean2((cs0 + cb0 + co20));
    Cell_SOC{1} = us_n;
    Cell_B{1} = ub_n;
    Cell_co2{1} = uco2_n;

    dec = ks .* us_n .* ub_n;
    dec_m(1)=mean2(dec);
    dec_std(1)=std2(dec);
    Cs_mean(1)=mean2(us_n);
    Cs_std(1)=std2(us_n);
    Cs_m4(1)=moment2(us_n,4);
    Cs_m3(1)=moment2(us_n,3);

    Cb_mean(1)=mean2(ub_n);
    Cb_std(1)=std2(ub_n);

    rho_CsCb(1)=corr2(us_n,ub_n);
    rho_Csks(1)=corr2(us_n,ks);
    rho_Cbks(1)=corr2(ub_n,ks);

    cov_cscb(1)= rho_CsCb(1)*std2(us_n)*std2(ub_n);
    cov_csks(1)= rho_Csks(1)*std2(us_n)*std2(ks);
    cov_cbks(1)= rho_Cbks(1)*std2(ub_n)*std2(ks);
    E_cs_cb_ks(1) = exp_CsCbks(us_n,ub_n,ks);
    
    total_cov(1) = (1-Y).*( cov_csks(1).*Cb_mean(1) + cov_cbks(1).*Cs_mean(1)...
                      + cov_cscb(1).*mean2(ks)+ E_cs_cb_ks(1));

    for n = 2:Nt
        [SOC, B, co2] = mult_diff_like_mass_Transfer_mex(us_n, ub_n, uco2_n, ks,...
            kb, Nx, Ny, Iteration, alpha(a1), dt, Imic, Y, omega, max_iter, tol, theta);


        % Update u_n before next step
        us_n = SOC;
        ub_n = B;
        uco2_n = co2;
        SOC = us_n;
        B = ub_n;
        co2 = uco2_n;
        sprintf('time step=%d', n)
        %==============post-proccessenig===================================
        mass_balance(n) = (sum(us_n(:)) + sum(ub_n(:)) + sum(uco2_n(:)) ...
            -((Imic) * dt * (n - 1) + (cs0 + cb0 + co20)) * Nx * Ny);

        Cell_SOC{n} = us_n;
        Cell_B{n} = ub_n;
        Cell_co2{n} = uco2_n;

        dec = ks .* us_n .* ub_n;
        dec_m(n)=mean2(dec);
        dec_std(n)=std2(dec);
        Cs_mean(n)=mean2(us_n);
        Cs_std(n)=std2(us_n);
        Cs_m4(n)=moment2(us_n,4);
        Cs_m3(n)=moment2(us_n,3);

        Cb_mean(n)=mean2(ub_n);
        Cb_std(n)=std2(ub_n);

        rho_CsCb(n)=corr2(us_n,ub_n);
        rho_Csks(n)=corr2(us_n,ks);
        rho_Cbks(n)=corr2(ub_n,ks);

        cov_cscb(n)= rho_CsCb(n)*std2(us_n)*std2(ub_n);
        cov_csks(n)= rho_Csks(n)*std2(us_n)*std2(ks);
        cov_cbks(n)= rho_Cbks(n)*std2(ub_n)*std2(ks);
        E_cs_cb_ks(n) = exp_CsCbks(us_n,ub_n,ks);

        total_cov(n) = (1-Y).*( cov_csks(n).*Cb_mean(n) + cov_cbks(n).*Cs_mean(n)...
                          + cov_cscb(n).*mean2(ks)+ E_cs_cb_ks(n));

    end
    save(['output\', 'Ph_mult_fullHetero_noCorr_',num2str(alpha(a1)), '.mat'])
end

toc
