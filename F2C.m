function [F2C] = F2C(solution,mesh)
% No top and bottom boundary condition since the u and v velocities here
% are equal to 0 and we just stay in the rows and shift across them
% Function file which creates the Divergence of Convective Flux matrix for
% the N-S equation that's a function of the 'u' velocity
u = solution.u;
v = solution.v;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
F2C = zeros(N,M);

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        %% Top
        F2C(N,j) = 0;
        
        % Bottom
        F2C(1,j) = 0;
        %% Main Domain
        F2C(i,j) = -((u(i,j+1)*v(i,j+1) - u(i,j-1)*v(i,j-1))/(2*mesh.dx)) - ((v(i+1,j)*v(i+1,j) - v(i-1,j)*v(i-1,j))/(2*mesh.dy));
        %% Left
        F2C(i,1) = -((u(i,2)*v(i,2) - u(i,M-1)*v(i,M-1))/(2*mesh.dx)) - ((v(i+1,1)*v(i+1,1) - v(i-1,1)*v(i-1,1))/(2*mesh.dy));  
        
        % Right
        F2C(i,M) = -((u(i,2)*v(i,2) - u(i,M-1)*v(i,M-1))/(2*mesh.dx)) - ((v(i+1,M)*v(i+1,M) - v(i-1,M)*v(i-1,M))/(2*mesh.dy));       
    end
end

% Top Left Corner
F2C(N,1) = 0;

% Bottom Left Corner
F2C(1,1) = 0;   

% Bottom Right Corner
F2C(1,M) = 0;

% Top Right Corner
F2C(N,M) = 0;
end