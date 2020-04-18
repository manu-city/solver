function [F1V1]=F1V1(solution,mesh,nonDimParams, F1V1)
% All bottom parts of the domain are 0 since there's no u velocity here
u = solution.u;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
Re = nonDimParams.reynolds_;

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Center
        F1V1(i,j) = (1/Re)*((u(i,j+1) - 2*u(i,j) + u(i,j-1))/((mesh.dx)^2));
        
        
        % Top
        F1V1(N,j) = (1/Re)*((u(N,j+1) - 2*u(N,j) + u(N,j-1))/((mesh.dx)^2));
    end 
end

for i = 2:N-1
    % Left
    F1V1(i,1) = (1/Re)*((u(i,2) - 2*u(i,1) + u(i,M-1))/((mesh.dx)^2));
    
    % Right
    F1V1(i,M) = (1/Re)*((u(i,2) - 2*u(i,M) + u(i,M-1))/((mesh.dx)^2));
end

% Top Left Corner
F1V1(N,1) = (1/Re)*((u(N,2) - 2*u(N,1) + u(N,M-1))/((mesh.dx)^2));

% Top Right Corner
F1V1(N,M) = (1/Re)*((u(N,2) - 2*u(N,M) + u(N,M-1))/((mesh.dx)^2));
end 