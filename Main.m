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
L                 = dimParams.o_y;
V                 = dimParams.bulkvel;

% Non-dimensionalisation of parameters
[nonDimParams]    = non_dim(dimParams, L, V);

%% Domain Discretisation 

% Number of x and y points
N                 = 100; 
M                 = 50; 

% Mesh the domain, obtain corresponding data from the mesh
[mesh]            = mesh(nonDimParams, N, M);

%% Stability

CFL = 0.6;
sigma = 0.3;

%% Initial Conditions

ic.u_velocity     = ones(N, M);  
ic.v_velocity     = ones(N, M);
ic.P              = ones(N, M);

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
    
    %% Divergence of Convective Flux matrix

    conv_flux.F1C = F1C(N,M,solution.u,solution.v,mesh);
    conv_flux.F2C = F2C(N,M,solution.u,solution.v,mesh);

    %% Divergence of Viscous Fluxes

    visc_flux.F1V1 = F1V1(N,M,solution.u,mesh,nonDimParams);
    visc_flux.F1V2 = F1V2(N,M,solution.u,mesh,nonDimParams);
    visc_flux.F2V1 = F2V1(N,M,solution.v,mesh,nonDimParams);
    visc_flux.F2V2 = F2V2(N,M,solution.v,mesh,nonDimParams);

    %% Divergence of Pressure Flux matrix

    press_flux.F1P = F1P(N,M,solution.P,mesh);
    press_flux.F2P = F2P(N,M,solution.P,mesh);

    %% G Functions

    G.G1 = conv_flux.F1C + visc_flux.F1V1 + force;
    G.G2 = conv_flux.F2C + visc_flux.F2V1;

    %% W Functions

    W.W1 = 0.5*(visc_flux.F1V2 + press_flux.F1P);
    W.W2 = 0.5*(visc_flux.F2V2 + press_flux.F2P);

    %% Prediction Velocity

    [pred] = prediction(mesh, nonDimParams, G, W, dt, t, k, solution, N, M);

    %% Inmersed Boundary Method

    [VOF, pred] = obstacles(mesh, nonDimParams, dimParams, N, M, pred);

    %% Poisson Solver

    [Pressure, A, b] = poisson(mesh, N, M, dt, solution);

    %% Update of velocity

    [solution] = update(mesh, Pressure, dt, pred, N, M);
    
    %% 
    
    UBulk = (sum(solution.u(2:end-1, 1)) + 0.5*solution.u(1,1) + 0.5*solution.u(N,1))/(N-1);
    duBulk = (1-UBulk)/dt;
    force = duBulk;
    
    %% Time marching

    [dt] = marching(mesh, solution, CFL, sigma, nonDimParams);

    k = k + 1;
    t(k) = t(k-1) + dt;
    
    %% Displaying
    
    contourf(flip(solution.u,1))
    colormap('jet')
    colorbar
    caxis([0 1])
    drawnow
end





