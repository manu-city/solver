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






