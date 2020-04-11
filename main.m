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
dimParams.clear   = (dimParams.L - dimParams.o_x * dimParams.o_num)/4;
dimParams.v       = 2e-5;

%% Domain Non-dimensionalisation 

% Length and velocity scales
L                 = dimParams.o_y;
V                 = dimParams.bulkvel;

% Non-dimensionalisation of parameters
[nonDimParams]    = non_dim(dimParams, L, V);

%% Domain Discretisation 

% Number of x and y points
N                 = 50;
M                 = 100;

% Mesh the domain, obtain corresponding data from the mesh
[mesh]            = mesh(nonDimParams, N, M);

%% N-S Discretisation

% Velocity vector
u                 = zeros(M,N);
v                 = zeros(M,N);

% Initial conditions; u(t=0) = bulkvel_, u(H_=0) = 0
u(1,:)            = nonDimParams.bulkvel_;
u(:,1)            = 0;

% Divergence of Convective Flux
[flux.F1C]        = F1C(N,M,u,v,mesh);
[flux.F2C]        = F2C(N,M,u,v,mesh);

% Divergence of Viscous Flux
[flux.F1V1]       = F1V1(N,M,u,mesh,nonDimParams);

[flux.F1V2]       = F1V2(N,M,u,mesh,nonDimParams);
[flux.F2V1]       = F2V1(N,M,v,mesh,nonDimParams);
[flux.F2V2]       = F2V2(N,M,v,mesh,nonDimParams);

% Divergence of Pressure Flux
P = zeros(M,N);
P(1,:) = 1;

[flux.F1P]        = F1P(N,M,P,mesh);
[flux.F2P]        = F2P(N,M,P,mesh);

% Summation of Fluxes
sum.G1            = flux.F1C + flux.F1V1;
sum.G2            = flux.F2C + flux.F2V2;

sum.W1            = 0.5*flux.F1V2 + 0.5*flux.F1P;
sum.W2            = 0.5*flux.F2V2 + 0.5*flux.F2P;