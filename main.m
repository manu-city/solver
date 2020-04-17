%%
clc
clear
close all
%% Geometry and domain features

% Dimensional parameters of the domain:
dimParams.L       = 1;
dimParams.H       = 0.25;
dimParams.o_x     = 0.05;
dimParams.o_y     = 0.05;
dimParams.o_num   = 4;
dimParams.bulkvel = 0.1;
dimParams.clear   = (dimParams.L - dimParams.o_x * dimParams.o_num)/dimParams.o_num;
dimParams.v       = 2e-5;

%% Domain Non-dimensionalisation 

% Length and velocity scales
L                 = dimParams.H;
V                 = dimParams.bulkvel;

% Non-dimensionalisation of parameters
[nonDimParams]    = non_dim(dimParams, L, V);

%% Domain Discretisation 

% Number of x and y points
N                 = 200;   % ROWS
M                 = 60;    % COLUMNS

% Mesh the domain, obtain corresponding data from the mesh
[mesh]            = mesh(nonDimParams, N, M);

%% Stability

CFL = 0.6;
sigma = 0.3;

%% Initial Conditions
N = mesh.ny;
M = mesh.nx;

ic.u_velocity     = ones(N, M);  
ic.v_velocity     = zeros(N, M);
ic.P              = zeros(N, M);

solution.u = ic.u_velocity;
solution.v = ic.v_velocity;
solution.P = ic.P;
force = 0;
%% Initial time step

t = zeros();
[dt] = marching(mesh, solution, CFL, sigma, nonDimParams);
k = 1;
t(k) = 0;

while true
    
    %% Divergence of Convective Flux matrix - 1
    conv_flux.F1C = F1C(solution,mesh);
    conv_flux.F2C = F2C(solution,mesh);
    %% Divergence of Viscous Fluxes - 2
    visc_flux.F1V1 = F1V1(solution,mesh,nonDimParams);
    visc_flux.F1V2 = F1V2(solution,mesh,nonDimParams);
    visc_flux.F2V1 = F2V1(solution,mesh,nonDimParams);
    visc_flux.F2V2 = F2V2(solution,mesh,nonDimParams);
    %% Divergence of Pressure Flux matrix - 3
    press_flux.F1P = F1P(solution,mesh);
    press_flux.F2P = F2P(solution,mesh);
    %% G Functions - 4
    G.G1 = conv_flux.F1C + visc_flux.F1V1 + force;
    G.G2 = conv_flux.F2C + visc_flux.F2V1;
    %% W Functions - 5
    W.W1 = 0.5*(visc_flux.F1V2 + press_flux.F1P);
    W.W2 = 0.5*(visc_flux.F2V2 + press_flux.F2P);
    %% Prediction Velocity - 6
    [pred] =  prediction(mesh, nonDimParams, G, W, dt, t, k, solution);
    %% Inmersed Boundary Method - 7
    [pred] = obstacles(mesh, nonDimParams, dimParams, pred);
    %% Poisson Solver - 8
    [Pressure, A, b] = poisson(mesh, dt, pred);
    %% Update of velocity - 9
    [solution] = update(mesh, Pressure, dt, pred);
    %% U Bulk - 10
    UBulk = (sum(solution.u(2:end-1, 1)) + 0.5*solution.u(1,1) + 0.5*solution.u(N,1))/(N);
    duBulk = (1-UBulk)/dt;
    force = duBulk;
    %% Time marching - 11
    [dt] = marching(mesh, solution, CFL, sigma, nonDimParams);

    k = k + 1;
    t(k) = t(k-1) + dt;
    %% Displaying - 12
    contourf((solution.u))
    title({['Iterations = ',num2str(k)]...
        [' Time Elapsed = ',num2str(t(end))]})
    colormap('parula(10)')
    colorbar
    caxis([min(solution.u,[],'all') max(solution.u,[],'all')])
    drawnow
    display(dt)
end