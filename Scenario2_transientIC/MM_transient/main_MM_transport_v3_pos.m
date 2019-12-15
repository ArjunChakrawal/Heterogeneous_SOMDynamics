
clc;
close all;
clearvars;

clc;close all;clearvars;
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


%% Simulate homogeneous Carbon dynamic for MM kinetics
load('Spatial_field\field_Csp.mat');
load('Spatial_field\field_Cb.mat');
cs0 = mean2(spCs);
cb0 = mean2(spCb);
co20 = 0;
Imic = cs0 / 200000; % external C input rate fgC/h
km = ((km / 1000) .* rho .* 1e-12 .* 1e15) .* Vol_mic; %fgC

c(1, :) = [cs0, cb0, co20];
St = ones(1, length(t)) .* Imic; % vector of external C input at each time step

f = @(tt, c)[-ksmm * c(1) * c(2) / (km + c(1)) + kb * c(2) + interp1(t, St, tt); ...
    Y * ksmm * c(1) * c(2) / (km + c(1)) - kb * c(2); ...
    (1 - Y) * ksmm * c(1) * c(2) / (km + c(1))];

options = odeset('Stats', 'off', 'AbsTol', 1e-9, 'RelTol', 1e-9); % MATLAB ode options
[tt, Ch] = ode45(f, t, c(1, :), options); %homogeneous solution

dec0 = ksmm .* Ch(:, 1) .* Ch(:, 2) ./ (km + Ch(:, 1)); %homogeneous rate of decomposition
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

%% Biophysical heterogeneity only i.e. only substrate and microbial C are spatially heterogeneous

% Relaxation coeff
omega = 0.1;
max_iter = 500;
tol = 1e-6;
abs_tol = 1e-12;
rel_tol = 1e-6;
% Crank–Nicolson theta parameter
theta = 0.5;

% hbar = parfor_progressbar(Nt,'time loop is running, please wait!');  %create the progress bar
mass_balance = zeros(Nt, length(alpha));

tic
for a1 = 1:length(alpha)
    %positive correlation
    ksmm = ones(Nx, Ny) .* ksmm;
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

    d2dec_dcs2 = @(csm, cbm, ksmmf, kmf)(-(2 .* cbm .* kmf .* ksmmf) ./ (csm + kmf).^3);
    d2dec_dkm2 = @(csm, cbm, ksmmf, kmf) (2 * cbm * csm * ksmmf) ./ (csm + kmf).^3;
    d2dec_dkmdksmm = @(csm, cbm, ksmmf, kmf) - (cbm * csm) ./ (csm + kmf).^2;
    d2dec_dcsdcb = @(csm, cbm, ksmmf, kmf)(ksmmf .* kmf) ./ ((csm + kmf).^2);
    d2dec_dcsdksmm = @(csm, cbm, ksmmf, kmf)(cbm * kmf) ./ (csm + kmf).^2;
    d2dec_dcbdksmm = @(csm, cbm, ksmmf, kmf)csm ./ (csm + kmf);
    d2dec_dcsdkm = @(csm, cbm, ksmmf, kmf) (cbm * ksmmf * (csm - kmf)) ./ (csm + kmf).^3;
    d2dec_dcbdkm = @(csm, cbm, ksmmf, kmf) - (csm * ksmmf) ./ (csm + kmf).^2;


    Cs_mean = zeros(length(t), 1);
    Cs_std = Cs_mean;
    Cs_m3 = Cs_mean;
    Cb_std = Cs_mean;
    CO2_mean = Cs_mean;
    CO2_std = Cs_mean;
    dec_m = Cs_mean;
    dec_std = Cs_mean;
    cov_cscb = Cs_mean;
    Cs_m4 = zeros(length(t), 1);
    Cb_mean = zeros(length(t), 1);
    rho_CsCb = zeros(length(t), 1);
    E_cs2_cb = zeros(length(t), 1);
    E_cs_cb2 = zeros(length(t), 1);
    E_cs2_cb2 = zeros(length(t), 1);

    rho_cs_ksmm = rho_CsCb;
    rho_cs_km = rho_CsCb;
    rho_cb_ksmm = rho_CsCb;
    rho_cb_km = rho_CsCb;
    rho_ksmm_km = rho_CsCb;

    cov_cs_ksmm = rho_cs_ksmm;
    cov_cs_km = rho_CsCb;
    cov_cb_ksmm = rho_CsCb;
    cov_cb_km = rho_CsCb;
    cov_ksmm_km = rho_CsCb;

%===========INITIAL CONDITION==========================================
    mass_balance(1) = mean2(us_n(:)) + mean2(ub_n(:)) + mean2(uco2_n(:)) ...
        -mean2((cs0 + cb0 + co20));
    Cell_SOC{1} = us_n;
    Cell_B{1} = ub_n;
    Cell_co2{1} = uco2_n;

    dec = ksmm .* us_n .* ub_n ./ (km + us_n);
    dec_m(1) = mean2(dec);
    dec_std(1) = std2(dec);

    Cs_mean(1) = mean2(us_n);
    Cs_std(1) = std2(us_n);
    Cs_m4(1) = moment2(us_n, 4);
    Cs_m3(1) = moment2(us_n, 3);

    Cb_mean(1) = mean2(ub_n);
    Cb_std(1) = std2(ub_n);

    CO2_mean(1) = mean2(uco2_n);
    CO2_std(1) = std2(uco2_n);

    rho_CsCb(1) = corr2(us_n, ub_n);
    rho_cs_ksmm(1) = corr2(us_n, ksmm);
    rho_cs_km(1) = corr2(us_n, km);
    rho_cb_ksmm(1) = corr2(ub_n, ksmm);
    rho_cb_km(1) = corr2(ub_n, km);
    rho_ksmm_km(1) = corr2(ksmm, km);


    cov_cscb(1) = rho_CsCb(1) * std2(us_n) * std2(ub_n);
    cov_cs_ksmm(1) = rho_cs_ksmm(1) * std2(us_n) * std2(ksmm);
    cov_cs_km(1) = rho_cs_km(1) * std2(us_n) * std2(km);
    cov_cb_ksmm(1) = rho_cb_ksmm(1) * std2(ub_n) * std2(ksmm);
    cov_cb_km(1) = rho_cb_km(1) * std2(ub_n) * std2(km);
    cov_ksmm_km(1) = rho_ksmm_km(1) * std2(ksmm) * std2(km);


    cov1(1) = 0.5 * d2dec_dcs2(Cs_mean(1), Cb_mean(1), mean2(ksmm), mean2(km)) * (std2(us_n))^2;
    cov2(1) = 0.5 * d2dec_dkm2(Cs_mean(1), Cb_mean(1), mean2(ksmm), mean2(km)) * (std2(km))^2;
    cov3(1) = d2dec_dkmdksmm(Cs_mean(1), Cb_mean(1), mean2(ksmm), mean2(km)) * cov_ksmm_km(1);

    cov4(1) = d2dec_dcsdcb(Cs_mean(1), Cb_mean(1), mean2(ksmm), mean2(km)) * cov_cscb(1);
    cov5(1) = d2dec_dcsdksmm(Cs_mean(1), Cb_mean(1), mean2(ksmm), mean2(km)) * cov_cs_ksmm(1);
    cov6(1) = d2dec_dcbdksmm(Cs_mean(1), Cb_mean(1), mean2(ksmm), mean2(km)) * cov_cb_ksmm(1);
    cov7(1) = d2dec_dcsdkm(Cs_mean(1), Cb_mean(1), mean2(ksmm), mean2(km)) * cov_cs_km(1);
    cov8(1) = d2dec_dcbdkm(Cs_mean(1), Cb_mean(1), mean2(ksmm), mean2(km)) * cov_cb_km(1);


    total_cov(1) = (1 - Y) .* (cov1(1) + cov2(1) + cov3(1) + cov4(1) + ... .
        cov5(1) + cov6(1) + cov7(1) + cov8(1));

    for n = 2:Nt
        [SOC, B, co2] = diff_like_mass_Transfer_mex(us_n, ub_n, uco2_n, ksmm, km, kb, Nx, Ny, ...
            Iteration, alpha(a1), dt, Imic, Y, omega, max_iter, tol, theta);

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

        dec = ksmm .* us_n .* ub_n ./ (km + us_n);
        dec_m(n) = mean2(dec);
        dec_std(n) = std2(dec);

        Cs_mean(n) = mean2(us_n);
        Cs_std(n) = std2(us_n);
        Cs_m4(n) = moment2(us_n, 4);
        Cs_m3(n) = moment2(us_n, 3);

        Cb_mean(n) = mean2(ub_n);
        Cb_std(n) = std2(ub_n);

        CO2_mean(n) = mean2(uco2_n);
        CO2_std(n) = std2(uco2_n);

        rho_CsCb(n) = corr2(us_n, ub_n);
        rho_cs_ksmm(n) = corr2(us_n, ksmm);
        rho_cs_km(n) = corr2(us_n, km);
        rho_cb_ksmm(n) = corr2(ub_n, ksmm);
        rho_cb_km(n) = corr2(ub_n, km);
        rho_ksmm_km(n) = corr2(ksmm, km);


        cov_cscb(n) = rho_CsCb(n) * std2(us_n) * std2(ub_n);
        cov_cs_ksmm(n) = rho_cs_ksmm(n) * std2(us_n) * std2(ksmm);
        cov_cs_km(n) = rho_cs_km(n) * std2(us_n) * std2(km);
        cov_cb_ksmm(n) = rho_cb_ksmm(n) * std2(ub_n) * std2(ksmm);
        cov_cb_km(n) = rho_cb_km(n) * std2(ub_n) * std2(km);
        cov_ksmm_km(n) = rho_ksmm_km(n) * std2(ksmm) * std2(km);


        cov1(n) = 0.5 * d2dec_dcs2(Cs_mean(n), Cb_mean(n), mean2(ksmm), mean2(km)) * (std2(us_n))^2;
        cov2(n) = 0.5 * d2dec_dkm2(Cs_mean(n), Cb_mean(n), mean2(ksmm), mean2(km)) * (std2(km))^2;
        cov3(n) = d2dec_dkmdksmm(Cs_mean(n), Cb_mean(n), mean2(ksmm), mean2(km)) * cov_ksmm_km(n);

        cov4(n) = d2dec_dcsdcb(Cs_mean(n), Cb_mean(n), mean2(ksmm), mean2(km)) * cov_cscb(n);
        cov5(n) = d2dec_dcsdksmm(Cs_mean(n), Cb_mean(n), mean2(ksmm), mean2(km)) * cov_cs_ksmm(n);
        cov6(n) = d2dec_dcbdksmm(Cs_mean(n), Cb_mean(n), mean2(ksmm), mean2(km)) * cov_cb_ksmm(n);
        cov7(n) = d2dec_dcsdkm(Cs_mean(n), Cb_mean(n), mean2(ksmm), mean2(km)) * cov_cs_km(n);
        cov8(n) = d2dec_dcbdkm(Cs_mean(n), Cb_mean(n), mean2(ksmm), mean2(km)) * cov_cb_km(n);


        total_cov(n) = (1 - Y) .* (cov1(n) + cov2(n) + cov3(n) + cov4(n) + ... .
            cov5(n) + cov6(n) + cov7(n) + cov8(n));
        E_cs2_cb(n) = exp_A2_B(us_n, ub_n);
        E_cs_cb2(n) = exp_A_B2(us_n, ub_n);
        E_cs2_cb2(n) = exp_A2_B2(us_n, ub_n);
    end
    save(['results\output_transport\varing_alpha\', 'results_v3_Ph_mm_pos_', ...
        num2str(alpha(a1)), '.mat'])

end

toc
%%
rmpath(selpath)
rmpath([selpath,'\Scenario2_transientIC'])
