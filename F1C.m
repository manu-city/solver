function [F1C] = F1C(solution,mesh, F1C)
% No bottom condition since all velocities at the bottom are equal to 0
% Function file which creates the Divergence of Convective Flux matrix for
% the N-S equation that's a function of the 'u' velocity
u = solution.u;
v = solution.v;
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)

for j = 2:M-1  
    % Top
    F1C(N,j) = -((u(N,j+1)*u(N,j+1) - u(N,j-1)*u(N,j-1))/(2*mesh.dx)) + ((u(N-1,j)*v(N-1,j))/(mesh.dy));% COLUMNS
    for i = 2:N-1     
        % Main Domain
        F1C(i,j) = -((u(i,j+1)*u(i,j+1) - u(i,j-1)*u(i,j-1))/(2*mesh.dx)) - ((u(i+1,j)*v(i+1,j) - u(i-1,j)*v(i-1,j))/(2*mesh.dy));
    end
end

for i = 2:N-1
    %% Left
    F1C(i,1) = -((u(i,2)*u(i,2) - -u(i,M-1)*-u(i,M-1))/(2*mesh.dx)) - ((u(i+1,1)*v(i+1,1) - u(i-1,1)*v(i-1,1))/(2*mesh.dy));
    
    % Right
    F1C(i,M) = -((u(i,2)*u(i,2) - u(i,M-1)*u(i,M-1))/(2*mesh.dx)) - ((u(i+1,M)*v(i+1,M) - u(i-1,M)*v(i-1,M))/(2*mesh.dy));
end

% Top Left Corner
F1C(N,1) = -((u(N,2)*u(N,2) - u(N,M-1)*u(N,M-1))/(2*mesh.dx)) + ((u(N-1,1)*v(N-1,1))/(mesh.dy));

% Top Right Corner
F1C(N,M) = -((u(N,2)*u(N,2) - u(N,M-1)*u(N,M-1))/(2*mesh.dx)) + ((u(N-1,M))*v(N-1,M))/(mesh.dy);
end       