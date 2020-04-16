function [F2V2]=F2V2(solution,mesh,nonDimParams)
v = solution.v;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
F2V2=zeros(N,M);
Re = nonDimParams.reynolds_;

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Center
        F2V2(i,j) = (1/Re)*((v(i+1,j) - 2*v(i,j) + v(i-1,j))/((mesh.dy)^2));
        
        % Left
        F2V2(i,1) = (1/Re)*((v(i+1,1) - 2*v(i,1) + v(i-1,1))/((mesh.dy)^2));
        
        % Right
        F2V2(i,M) = (1/Re)*((v(i+1,M) - 2*v(i,M) + v(i-1,M))/((mesh.dy)^2));
        
        % Bottom
        F2V2(N,j) = 0;
        
        % Top
        F2V2(1,j) = 0;
    end 
end 
        % Bottom Left Corner
        F2V2(N,1) = 0;
        
        % Bottom Right Corner
        F2V2(N,M) = 0;
        
        % Top Left Corner
        F2V2(1,1) = 0;
        
        % Top Right Corner
        F2V2(1,M) = 0;
end 