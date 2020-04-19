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
N                 = 100;   % ROWS
M                 = 60;    % COLUMNS

% Mesh the domain, obtain corresponding data from the mesh
[mesh]            = mesh(nonDimParams, N, M);

%% Stability

CFL = 0.6;
sigma = 0.3;

%% Initial Conditions
% Used the mesh values since we should be using gridpoints instead of the
% cells
% Changed the orientatin of the matrix so that now the bottom of the domain
% will be at the top of the matrix and the top of the domain at the bottom
% of the matrix, this is constant throughout now
N = mesh.ny;
M = mesh.nx;

ic.u_velocity     = ones(N, M);
ic.v_velocity     = zeros(N, M);
ic.P              = zeros(N, M);

solution.u        = ic.u_velocity;
solution.v        = ic.v_velocity;
solution.P        = ic.P;
force             = 0;

%% Initialising space

% Fluxes
conv_flux.F1C     = zeros(N,M);
conv_flux.F2C     = zeros(N,M);
visc_flux.F1V1    = zeros(N,M);
visc_flux.F1V2    = zeros(N,M);
visc_flux.F2V1    = zeros(N,M);
visc_flux.F2V2    = zeros(N,M);
press_flux.F1P    = zeros(N,M);
press_flux.F2P    = zeros(N,M);

% Poisson Matrix and vector
A                 = zeros(N*M);
b                 = zeros(N*M,1);
[AP]              = Apoisson(mesh, A);

% u and v star and corresponding A matrices
pred.predU        = zeros(N,M);
pred.predV        = zeros(N,M);
predMat.A1        = zeros(N-1,1);
predMat.A1_       = zeros(N-2,1);

% Vortices and Streamlines
dS                = 5;  % Spacing between streamlines
dV                = 2;  % Spacing between vortices
[sx,sy]           = meshgrid(0:dS:M-1,0:dS:N-1);
[vx,vy]           = meshgrid(0:dV:M-1,0:dV:N-1);
%% Initial time step

t = zeros();
[dt] = marching(mesh, solution, CFL, sigma, nonDimParams);
k = 1;
t(k) = 0;

while true
    
    %% Divergence of Convective Flux matrix - 1
    
    conv_flux.F1C = F1C(solution,mesh, conv_flux.F1C);
    conv_flux.F2C = F2C(solution,mesh, conv_flux.F2C);
    
    %% Divergence of Viscous Fluxes - 2
    
    visc_flux.F1V1 = F1V1(solution,mesh,nonDimParams,visc_flux.F1V1);
    visc_flux.F1V2 = F1V2(solution,mesh,nonDimParams,visc_flux.F1V2);
    visc_flux.F2V1 = F2V1(solution,mesh,nonDimParams,visc_flux.F2V1);
    visc_flux.F2V2 = F2V2(solution,mesh,nonDimParams,visc_flux.F2V2);
    
    %% Divergence of Pressure Flux matrix - 3
    
    press_flux.F1P = F1P(solution,mesh,press_flux.F1P);
    press_flux.F2P = F2P(solution,mesh,press_flux.F2P);
    
    %% G Functions - 4
    
    G.G1 = conv_flux.F1C + visc_flux.F1V1 + force;
    G.G2 = conv_flux.F2C + visc_flux.F2V1;
    
    %% W Functions - 5
    
    W.W1 = 0.5*(visc_flux.F1V2 + press_flux.F1P);
    W.W2 = 0.5*(visc_flux.F2V2 + press_flux.F2P);
    
    %% Prediction Velocity - 6
    
    [pred] =  prediction(mesh, nonDimParams, G, W, dt, t, k, solution, pred, predMat);
    
    %% Inmersed Boundary Method - 7
    
    [pred] = obstacles(mesh, nonDimParams, dimParams, pred);
   
    %% Poisson Solver - 8
    
    [Pressure] = poisson(mesh, dt, pred, AP, b);
    
    %% Update of velocity - 9
    
    [solution] = update(mesh, Pressure, dt, pred, solution);
    
    %% U Bulk - 10
    
    UBulk = (sum(solution.u(2:end-1, 1)) + 0.5*solution.u(1,1) + 0.5*solution.u(N,1))/(N);
    duBulk = (1-UBulk)/dt;
    force = duBulk;
    
    %% Time marching - 11
    
    [dt] = marching(mesh, solution, CFL, sigma, nonDimParams);

    k = k + 1;
    t(k) = t(k-1) + dt;
    
    %% Displaying - Contours
    figure(1)
    contourf((solution.u))
    title('X-Componenet [U]')
    colormap('parula(10)')
    colorbar
    caxis([min(solution.u,[],'all') max(solution.u,[],'all')])
    
    figure(2)
    contourf((solution.v))
    title('Y-Componenet [V]')
    colormap('parula(10)')
    colorbar
    caxis([min(solution.v,[],'all') max(solution.v,[],'all')])
    
    drawnow
    %% Displaying - Vorticies
    figure(3)
    hold on
    cav = curl(vx,vy,solution.u(1:dV:end,1:dV:end),solution.v(1:dV:end,1:dV:end));
    pcolor(vx,vy,cav)
    shading interp
    hold on
    quiver(vx,vy,solution.u(1:dV:end,1:dV:end),solution.v(1:dV:end,1:dV:end),'k')
    title('Vorticies')
    xlim manual
    ylim manual
    axis([0 M-1 0 N-1])
    colormap('parula(10)')
    colorbar
    hold off
    
    drawnow
    %% Displaying - Streamlines
    figure(4)
    streamline(solution.u,solution.v,sx,sy);
    title('Streamlines')
    xlim manual
    ylim manual
    axis([dS M-1 dS N-1])
    hold off
    
    drawnow
    %% Displaying - Iteration & Time Step
    
    fprintf('Iterations = %1.0f and Time Step = %4.5f\n',k-1,t(end));
    
end