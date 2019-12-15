clc;
close all;
clearvars;

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

% fraction of rate of susbtrate transferred to neigbhour grid cells
alpha = 0;
%% Estimation of kinetic parameter

load('C:\Users\arch9809\Box Sync\Stockholm Unviersity\Notebook\Papaer1\gmd\review\code\Scenario2_transientIC\Spatial_field\field_Csp.mat');
load('C:\Users\arch9809\Box Sync\Stockholm Unviersity\Notebook\Papaer1\gmd\review\code\Scenario2_transientIC\Spatial_field\field_Cb.mat');
cs0 = mean2(spCs);cb0 = mean2(spCb);
Imic = cs0 / 200000; % external C input rate fgC/h
Km_inv=cb0*1e-12./(rho*Vmic_cc)*8;
ksmm_inv=ksmm/4;
Km_inv = ((Km_inv / 1000) .* rho .* 1e-12 .* 1e15) .* Vol_mic; %fgC
%%
load('C:\Users\arch9809\Box Sync\Stockholm Unviersity\Notebook\Papaer1\gmd\review\code\MM_inv_Steady_State\Spatial_field\ph_neg_cs_mm_inv.mat');
load('C:\Users\arch9809\Box Sync\Stockholm Unviersity\Notebook\Papaer1\gmd\review\code\MM_inv_Steady_State\Spatial_field\ph_cb.mat');

% Relaxation coeff
omega = 0.1;
max_iter = 500;
tol = 1e-9;
abs_tol = 1e-12;
rel_tol = 1e-6;
% Crank–Nicolson theta parameter
theta = 0.5;

%positive correlation
ksmm_inv = ones(Nx, Ny) .* ksmm_inv;
kb = ones(Nx, Ny) .* kb;
Km_inv = ones(Nx, Ny) .* Km_inv;
% The initial condition is easy to fill in:
us_n = spCs;
ub_n = spCb;
uco2_n = zeros(Nx, Ny);


Iteration='SOR';
tic
[SOC, B, co2, totalC] = MM_inv_diff_like_mass_Transfer(us_n, ub_n, uco2_n, ksmm_inv, Km_inv,...
    kb, Nx, Ny, Iteration, alpha, dt, Imic, Y, omega, max_iter, tol, theta);
toc

tic
[SOC, B, co2, totc] = MM_inv_diff_like_mass_Transfer_mex(us_n, ub_n, uco2_n, ksmm_inv, Km_inv,...
    kb, Nx, Ny, Iteration, alpha, dt, Imic, Y, omega, max_iter, tol, theta);

toc

