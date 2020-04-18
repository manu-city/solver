function [F2V1]=F2V1(solution,mesh,nonDimParams, F2V1)
% No v velocity in top or bottom of domain so these boundaries all equal to
% zero, we're going across these rows so there would be no change from zero
% at the top and bottom boundary
v = solution.v;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
Re = nonDimParams.reynolds_;

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Center
        F2V1(i,j) = (1/Re)*((v(i,j+1) - 2*v(i,j) + v(i,j-1))/((mesh.dx)^2));        
    end  
end 
for i = 2:N-1        
    % Left
    F2V1(i,1) = (1/Re)*((v(i,2) - 2*v(i,1) + v(i,M-1))/((mesh.dx)^2));
    
    % Right
    F2V1(i,M) = (1/Re)*((v(i,2) - 2*v(i,M) + v(i,M-1))/((mesh.dx)^2));
end

end 