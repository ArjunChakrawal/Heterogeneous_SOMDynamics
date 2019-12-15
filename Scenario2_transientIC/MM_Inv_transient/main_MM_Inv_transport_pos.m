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
%%
load('Spatial_field\field_Csp.mat');
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
[~, Ch2] = ode45(f, t, c(1, :), options);

Cd1 = Ch2 .* 1e-12 ./ (rho * Vmic_cc); %Time evolution of carbon mg C at microbsite

decmminv = ksmm_inv .* Ch2(:, 1) .* Ch2(:, 2) ./ (Km_inv + Ch2(:, 2));

figure(1);
p2 = plot(t./24, Cd1(:, 1), '--');
hold on;
plot(t./24, Cd1(:, 2), '--');
plot(t./24, Cd1(:, 3), '--');
xlabel('Time(day)');
ylabel('(mg C/gSoil) ');
legend('Substrate', 'Biomass', 'CO_2'); title('Homogeneous case')
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', ...
    'LineWidth', 0.5, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02, .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XColor', [1, 1, 1]*0.3, 'YColor', [1, 1, 1]*0.3);
set(gcf, 'Color', 'w')

figure; plot(t./24, (1 - Y).*decmminv.*1e-12./(rho * Vmic_cc))
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
rel_tol = 1e-12;
% Crank–Nicolson theta parameter
theta = 0.5;

% hbar = parfor_progressbar(Nt,'time loop is running, please wait!');  %create the progress bar
mass_balance = zeros(Nt, length(alpha));

tic
for a1 = 1:length(alpha)
    %positive correlation
    ksmm_inv = ones(Nx, Ny) .* ksmm_inv;
    kb = ones(Nx, Ny) .* kb;
    Km_inv = ones(Nx, Ny) .* Km_inv;
    % The initial condition is easy to fill in:
    us_n = spCs;
    ub_n = spCb;
    uco2_n = zeros(Nx, Ny);

    SOC = zeros(Nx, Ny); % unknown SOC at new time level
    B = zeros(Nx, Ny);
    co2 = zeros(Nx, Ny);

    Cell_SOC = cell(Nt, 1);
    Cell_B = Cell_SOC;
    Cell_co2 = Cell_SOC;
    totalCarbon = zeros(Nt, 1);
    d2dec_dcb2 = @(csm, cbm, ksmmf, kmf)(-(2 .* csm .* kmf .* ksmmf) ./ (cbm + kmf).^3);
    d2dec_dKm_inv2 = @(csm, cbm, ksmmf, kmf) (2 * cbm * csm * ksmmf) ./ (cbm + kmf).^3;
    d2dec_dKm_invdksmm_inv = @(csm, cbm, ksmmf, kmf) - (cbm * csm) ./ (cbm + kmf).^2;
    d2dec_dcsdcb = @(csm, cbm, ksmmf, kmf)(ksmmf .* kmf) ./ ((cbm + kmf).^2);
    d2dec_dcsdksmm_inv = @(csm, cbm, ksmmf, kmf) cbm ./ (cbm + kmf);
    d2dec_dcbdksmm_inv = @(csm, cbm, ksmmf, kmf)(csm .* kmf) ./ (cbm + kmf).^2;
    d2dec_dcsdKm_inv = @(csm, cbm, ksmmf, kmf) - (cbm * ksmmf) ./ (cbm + kmf).^2;
    d2dec_dcbdKm_inv = @(csm, cbm, ksmmf, kmf) (csm * ksmmf * (cbm - kmf)) ./ (cbm + kmf).^3;


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

    rho_cb_ksmm_inv = rho_CsCb;
    rho_cs_Km_inv = rho_CsCb;
    rho_cs_ksmm_inv = rho_CsCb;
    rho_cb_Km_inv = rho_CsCb;
    rho_ksmm_inv_Km_inv = rho_CsCb;

    cov_cs_ksmm_inv = rho_CsCb;
    cov_cs_Km_inv = rho_CsCb;
    cov_cb_ksmm_inv = rho_CsCb;
    cov_cb_Km_inv = rho_CsCb;
    cov_ksmm_inv_Km_inv = rho_CsCb;
    

    %===========INITIAL CONDITION==========================================
    mass_balance(1, a1) = mean2(us_n(:)) + mean2(ub_n(:)) + mean2(uco2_n(:)) ...
        -mean2((cs0 + cb0 + co20));
    totalCarbon(1) = mean2((cs0 + cb0 + co20))*Nx*Ny;
    Cell_SOC{1} = us_n;
    Cell_B{1} = ub_n;
    Cell_co2{1} = uco2_n;
    dec = ksmm_inv .* us_n .* ub_n ./ (Km_inv + ub_n);


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
    rho_cs_ksmm_inv(1) = corr2(us_n, ksmm_inv);
    rho_cs_Km_inv(1) = corr2(us_n, Km_inv);
    rho_cb_ksmm_inv(1) = corr2(ub_n, ksmm_inv);
    rho_cb_Km_inv(1) = corr2(ub_n, Km_inv);
    rho_ksmm_inv_Km_inv(1) = corr2(ksmm_inv, Km_inv);


    cov_cscb(1) = rho_CsCb(1) * std2(us_n) * std2(ub_n);
    cov_cs_ksmm_inv(1) = rho_cs_ksmm_inv(1) * std2(us_n) * std2(ksmm_inv);
    cov_cs_Km_inv(1) = rho_cs_Km_inv(1) * std2(us_n) * std2(Km_inv);
    cov_cb_ksmm_inv(1) = rho_cb_ksmm_inv(1) * std2(ub_n) * std2(ksmm_inv);
    cov_cb_Km_inv(1) = rho_cb_Km_inv(1) * std2(ub_n) * std2(Km_inv);
    cov_ksmm_inv_Km_inv(1) = rho_ksmm_inv_Km_inv(1) * std2(ksmm_inv) * std2(Km_inv);


    cov1(1) = 0.5 * d2dec_dcb2(Cs_mean(1), Cb_mean(1), mean2(ksmm_inv), mean2(Km_inv)) * (std2(ub_n))^2;
    cov2(1) = 0.5 * d2dec_dKm_inv2(Cs_mean(1), Cb_mean(1), mean2(ksmm_inv), mean2(Km_inv)) * (std2(Km_inv))^2;
    cov3(1) = d2dec_dKm_invdksmm_inv(Cs_mean(1), Cb_mean(1), mean2(ksmm_inv), mean2(Km_inv)) * cov_ksmm_inv_Km_inv(1);

    cov4(1) = d2dec_dcsdcb(Cs_mean(1), Cb_mean(1), mean2(ksmm_inv), mean2(Km_inv)) * cov_cscb(1);
    cov5(1) = d2dec_dcsdksmm_inv(Cs_mean(1), Cb_mean(1), mean2(ksmm_inv), mean2(Km_inv)) * cov_cs_ksmm_inv(1);
    cov6(1) = d2dec_dcbdksmm_inv(Cs_mean(1), Cb_mean(1), mean2(ksmm_inv), mean2(Km_inv)) * cov_cb_ksmm_inv(1);
    cov7(1) = d2dec_dcsdKm_inv(Cs_mean(1), Cb_mean(1), mean2(ksmm_inv), mean2(Km_inv)) * cov_cs_Km_inv(1);
    cov8(1) = d2dec_dcbdKm_inv(Cs_mean(1), Cb_mean(1), mean2(ksmm_inv), mean2(Km_inv)) * cov_cb_Km_inv(1);


    total_cov(1) = (1 - Y) .* (cov1(1) + cov2(1) + cov3(1) + cov4(1) + ... .
        cov5(1) + cov6(1) + cov7(1) + cov8(1));


    for n = 2:Nt

        [SOC, B, co2, totc] = MM_inv_diff_like_mass_Transfer_mex(us_n, ub_n, uco2_n, ksmm_inv, Km_inv, ...
            kb, Nx, Ny, Iteration, alpha(a1), dt, Imic, Y, omega, max_iter, tol, theta);

        % Update u_n before next step
        us_n = SOC;
        ub_n = B;
        uco2_n = co2;
        SOC = us_n;
        B = ub_n;
        co2 = uco2_n;
        %==============post-proccessenig===================================
        mass_balance(n, a1) = sum(SOC(:)) + sum(B(:)) + sum(co2(:)) ...
            -Imic * dt * (n-1) * Nx * Ny - sum((cs0(:) + cb0(:))) * Nx * Ny;
        totalCarbon(n) = totc;

        Cell_SOC{n} = SOC;
        Cell_B{n} = B;
        Cell_co2{n} = co2;
        dec = ksmm_inv .* SOC .* B ./ (Km_inv + B);


        dec_m(n) = mean2(dec);
        dec_std(n) = std2(dec);

        Cs_mean(n) = mean2(SOC);
        Cs_std(n) = std2(SOC);
        Cs_m4(n) = moment2(SOC, 4);
        Cs_m3(n) = moment2(SOC, 3);

        Cb_mean(n) = mean2(B);
        Cb_std(n) = std2(B);

        CO2_mean(n) = mean2(co2);
        CO2_std(n) = std2(co2);

        rho_CsCb(n) = corr2(SOC, B);
        rho_cs_ksmm_inv(n) = corr2(SOC, ksmm_inv);
        rho_cs_Km_inv(n) = corr2(SOC, Km_inv);
        rho_cb_ksmm_inv(n) = corr2(B, ksmm_inv);
        rho_cb_Km_inv(n) = corr2(B, Km_inv);
        rho_ksmm_inv_Km_inv(n) = corr2(ksmm_inv, Km_inv);


        cov_cscb(n) = rho_CsCb(n) * std2(SOC) * std2(B);
        cov_cs_ksmm_inv(n) = rho_cs_ksmm_inv(n) * std2(SOC) * std2(ksmm_inv);
        cov_cs_Km_inv(n) = rho_cs_Km_inv(n) * std2(SOC) * std2(Km_inv);
        cov_cb_ksmm_inv(n) = rho_cb_ksmm_inv(n) * std2(B) * std2(ksmm_inv);
        cov_cb_Km_inv(n) = rho_cb_Km_inv(n) * std2(B) * std2(Km_inv);
        cov_ksmm_inv_Km_inv(n) = rho_ksmm_inv_Km_inv(n) * std2(ksmm_inv) * std2(Km_inv);


        cov1(n) = 0.5 * d2dec_dcb2(Cs_mean(n), Cb_mean(n), mean2(ksmm_inv), mean2(Km_inv)) * (std2(B))^2;
        cov2(n) = 0.5 * d2dec_dKm_inv2(Cs_mean(n), Cb_mean(n), mean2(ksmm_inv), mean2(Km_inv)) * (std2(Km_inv))^2;
        cov3(n) = d2dec_dKm_invdksmm_inv(Cs_mean(n), Cb_mean(n), mean2(ksmm_inv), mean2(Km_inv)) * cov_ksmm_inv_Km_inv(n);

        cov4(n) = d2dec_dcsdcb(Cs_mean(n), Cb_mean(n), mean2(ksmm_inv), mean2(Km_inv)) * cov_cscb(n);
        cov5(n) = d2dec_dcsdksmm_inv(Cs_mean(n), Cb_mean(n), mean2(ksmm_inv), mean2(Km_inv)) * cov_cs_ksmm_inv(n);
        cov6(n) = d2dec_dcbdksmm_inv(Cs_mean(n), Cb_mean(n), mean2(ksmm_inv), mean2(Km_inv)) * cov_cb_ksmm_inv(n);
        cov7(n) = d2dec_dcsdKm_inv(Cs_mean(n), Cb_mean(n), mean2(ksmm_inv), mean2(Km_inv)) * cov_cs_Km_inv(n);
        cov8(n) = d2dec_dcbdKm_inv(Cs_mean(n), Cb_mean(n), mean2(ksmm_inv), mean2(Km_inv)) * cov_cb_Km_inv(n);


        total_cov(n) = (1 - Y) .* (cov1(n) + cov2(n) + cov3(n) + cov4(n) + ... .
            cov5(n) + cov6(n) + cov7(n) + cov8(n));


        sprintf('time step=%d', n)
    end
    save(['output\', 'Ph_mm_inv_pos_', num2str(alpha(a1)), '.mat'])
end

toc



