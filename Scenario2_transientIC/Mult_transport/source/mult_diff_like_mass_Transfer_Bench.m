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
Nt = round(T/dt);
t = linspace(0, Nt*dt, Nt); % Mesh points in time

% fraction of rate of susbtrate transferred to neigbhour grid cells
alpha = [0, 0.01, 0.05, 0.1, 0.2];

%% %% ksmult estimation and Simulate homogeneous Carbon dynamic for MM kinetics
load('C:\Users\arch9809\Box Sync\Stockholm Unviersity\Notebook\Papaer1\gmd\review\code\Scenario2_transientIC\Spatial_field\field_Csn.mat');
load('C:\Users\arch9809\Box Sync\Stockholm Unviersity\Notebook\Papaer1\gmd\review\code\Scenario2_transientIC\Spatial_field\field_Cb.mat');
cs0=mean2(spCs)*1e-12./(rho*Vmic_cc);cb0=mean2(spCb)*1e-12./(rho*Vmic_cc);co20=0;
Imic=cs0/200000; T=24*1000;
ks=ks_mult_estimation(T,cs0,cb0,co20,kb,ksmm, km, Imic,Y,rho,Vol_mic);
%% Relaxation coeff
omega = 0.1;
max_iter = 500;
tol = 1e-6;
abs_tol = 1e-12;
rel_tol = 1e-6;
% Crank–Nicolson theta parameter
theta = 0.5;


a1 = 1;
%positive correlation
ks = ones(Nx, Ny) .* ks;
kb = ones(Nx, Ny) .* kb;
km = ones(Nx, Ny) .* km;

% The initial condition is easy to fill in:
us_n = spCs;
ub_n = spCb;
uco2_n = zeros(Nx, Ny);

Iteration='SOR';
tic
[SOC, B, co2] = mult_diff_like_mass_Transfer(us_n, ub_n, uco2_n, ks,...
            kb, Nx, Ny, Iteration, alpha(1), dt, Imic, Y, omega, max_iter, tol, theta);
toc

tic
[SOC, B, co2] =  mult_diff_like_mass_Transfer_mex(us_n, ub_n, uco2_n, ks,...
            kb, Nx, Ny, Iteration, alpha(1), dt, Imic, Y, omega, max_iter, tol, theta);
toc

