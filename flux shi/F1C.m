function [F1C] = F1C(N,M,u,v,mesh)
% Function file which creates the Divergence of Convective Flux matrix for
% the N-S equation that's a function of the 'u' velocity
u = solution.u;
v = solution.v;
F1C = zeros(N,M);

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        %% Top
        F1C(1,j) = -((u(1,j+1)*u(1,j+1) - u(1,j-1)*u(1,j-1))/(2*mesh.dx)) - ((u(2,j)*v(2,j))/(mesh.dy));
         
        % Bottom
        F1C(N,j) = - ((-u(N-1,j)*v(N-1,j))/(mesh.dy));
        %% Main Domain
        F1C(i,j) = -((u(i,j+1)*u(i,j+1) - u(i,j-1)*u(i,j-1))/(2*mesh.dx)) - ((u(i+1,j)*v(i+1,j) - u(i-1,j)*v(i-1,j))/(2*mesh.dy));
        %% Left
        F1C(i,1) = -((u(i,2)*u(i,2)) - (-u(i,M-1)*-u(i,M-1)))/(2*mesh.dx) -((u(i+1,1)*v(i+1,1)) - (u(i-1,1)*v(i-1,1)))/(2*mesh.dy); 
        
        % Right
        F1C(i,M) = -((u(i,2)*u(i,2) - u(i,M-1)*u(i,M-1))/(2*mesh.dx)) - ((u(i+1,M)*v(i+1,M) - u(i-1,M)*v(i-1,M))/(2*mesh.dy));       
    end
end

% Top Left Corner
F1C(1,1) = -((u(1,2)*u(1,2) - u(1,M-1)*u(1,M-1))/(2*mesh.dx)) - ((u(2,1)*v(2,1))/(mesh.dy));

% Bottom Left Corner
F1C(N,1) = - ((-u(N-1,1)*v(N-1,1))/(mesh.dy));

% Bottom Right Corner
F1C(N,M) = - ((-u(N-1,M)*v(N-1,M))/(mesh.dy));

% Top Right Corner
F1C(1,M) = -((u(1,2)*u(N,2) - u(1,M-1)*u(1,M-1))/(2*mesh.dx)) - ((u(2,M))/(mesh.dy));
end       