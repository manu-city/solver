function [F1V1]=F1V1(solution,mesh,nonDimParams)
u = solution.u;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
F1V1=zeros(N,M);
Re = nonDimParams.reynolds_;

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Center
        F1V1(i,j) = (1/Re)*((u(i,j+1) - 2*u(i,j) + u(i,j-1))/((mesh.dx)^2));
        
        % Left
        F1V1(i,1) = (1/Re)*((u(i,2) - 2*u(i,1) + u(i,M-1))/((mesh.dx)^2));
        
        % Right
        F1V1(i,M) = (1/Re)*((u(i,2) - 2*u(i,M) + u(i,M-1))/((mesh.dx)^2));
        
        % Bottom
        F1V1(1,j) = 0;
        
        % Top
        F1V1(N,j) = (1/Re)*((u(N,j+1) - 2*u(N,j) + u(N,j-1))/((mesh.dx)^2));
    end 
end 
        % Bottom Left Corner
        F1V1(1,1) = 0;
        
        % Bottom Right Corner
        F1V1(1,M) = 0;
        
        % Top Left Corner
        F1V1(N,1) = (1/Re)*((u(N,2) - 2*u(N,1) + u(N,M-1))/((mesh.dx)^2));
        
        % Top Right Corner
        F1V1(N,M) = (1/Re)*((u(N,2) - 2*u(N,M) + u(N,M-1))/((mesh.dx)^2));
end 