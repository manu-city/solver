function [F2V1]=F2V1(solution,mesh,nonDimParams)
v = solution.v;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
F2V1=zeros(N,M);
Re = nonDimParams.reynolds_;

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Center
        F2V1(i,j) = (1/Re)*((v(i,j+1) - 2*v(i,j) + v(i,j-1))/((mesh.dx)^2));
        
        % Left
        F2V1(i,1) = (1/Re)*((v(i,2) - 2*v(i,1) + v(i,M-1))/((mesh.dx)^2));
        
        % Right
        F2V1(i,M) = (1/Re)*((v(i,2) - 2*v(i,M) + v(i,M-1))/((mesh.dx)^2));
        
        % Bottom
        F2V1(N,j) = 0;
        
        % Top
        F2V1(1,j) = 0;
    end 
end 
        % Bottom Left Corner
        F2V1(N,1) = 0;
        
        % Bottom Right Corner
        F2V1(N,M) = 0;
        
        % Top Left Corner
        F2V1(1,1) = 0;
        
        % Top Right Corner
        F2V1(1,M) = 0;
end 