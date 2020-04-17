function [mesh] = mesh(nonDimParams, N, M)

%% 

% Can't figure out a way of doing all of this better.

%%

% dx and dy for uniform grid
mesh.dx = nonDimParams.L_/M;
mesh.dy = nonDimParams.H_/N;

% Number of index points
mesh.nx = nonDimParams.L_/mesh.dx + 1;
mesh.ny = nonDimParams.H_/mesh.dy + 1;

mesh.Dx = zeros();
mesh.Dy = zeros();

% vector containing all dx and dy positions
for i = 1:mesh.nx
    mesh.Dx(i, 1) = (i-1) * mesh.dx;
end

for i = 1:mesh.ny
    mesh.Dy(i, 1) = (i-1) * mesh.dy;
end

