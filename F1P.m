function [F1P] = F1P(solution,mesh, F1P)
P = solution.P;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
% Function file which creates the Divergence of Pressure Flux matrix for
% the N-S equation that's a function of the 'u' velocity


for j = 2:M-1         
    % Bottom
    F1P(1,j) = -(P(1,j+1) - P(1,j-1))/(2*mesh.dx);
    for i = 2:N-1     
        % Main
        F1P(i,j) = -(P(i,j+1) - P(i,j-1))/(2*mesh.dx);
        
        % Top
        F1P(N,j) = -(P(N,j+1) - P(N,j-1))/(2*mesh.dx);
    end
end

for i = 2:N-1
    % Left
    F1P(i,1) = -(P(i,2) - P(i,M-1))/(2*mesh.dx);
    
    % Right
    F1P(i,M) = -(P(i,2) - P(i,M-1))/(2*mesh.dx);
end
% Bottom Left
F1P(1,1) = -(P(1,2) - P(1,M-1))/(2*mesh.dx);

% Bottom Right
F1P(1,M) = -(P(1,2) - P(1,M-1))/(2*mesh.dx);

% Top Left
F1P(N,1) = -(P(N,2) - P(N,M-1))/(2*mesh.dx);

% Top Right
F1P(N,M) = -(P(N,2) - P(N,M-1))/(2*mesh.dx);
end