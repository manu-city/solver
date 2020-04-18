function [F1V2]=F1V2(solution,mesh,nonDimParams, F1V2)
% Similar to F1V1, bottom is all 0 since there's no u velocity
u = solution.u;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)

Re = nonDimParams.reynolds_;

for j = 2:M-1         
    % Top
    F1V2(N,j) = (1/Re)*((2*u(N-1,j) - 2*u(N,j))/((mesh.dy)^2));
    for i = 2:N-1     
        % Main
        F1V2(i,j) = (1/Re)*((u(i+1,j) - 2*u(i,j) + u(i-1,j))/((mesh.dy)^2));
    end 
end 

for i = 2:N-1
    % Left
    F1V2(i,1) = (1/Re)*((u(i+1,1) - 2*u(i,1) + u(i-1,1))/((mesh.dy)^2));
    
    % Right
    F1V2(i,M) = (1/Re)*((u(i+1,M) - 2*u(i,M) + u(i-1,M))/((mesh.dy)^2)); 
end

% Top Left Corner
F1V2(N,1) = (1/Re)*((2*u(N-1,j) - 2*u(N,1))/((mesh.dy)^2));

% Top Right Corner
F1V2(N,M) = (1/Re)*((2*u(N-1,j) - 2*u(N,M))/((mesh.dy)^2));
end 