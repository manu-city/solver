function [F2V2]=F2V2(solution,mesh,nonDimParams, F2V2)
% v is 0 at the top and bottom boundary
v = solution.v;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
Re = nonDimParams.reynolds_;

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Center
        F2V2(i,j) = (1/Re)*((v(i+1,j) - 2*v(i,j) + v(i-1,j))/((mesh.dy)^2));
    end 
end 

for i = 2:N-1
    % Left
    F2V2(i,1) = (1/Re)*((v(i+1,1) - 2*v(i,1) + v(i-1,1))/((mesh.dy)^2));
    
    % Right
    F2V2(i,M) = (1/Re)*((v(i+1,M) - 2*v(i,M) + v(i-1,M))/((mesh.dy)^2));
end
end 