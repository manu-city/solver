function [F1C] = F1C(N,M,u,v,mesh)
% Function file which creates the Divergence of Convective Flux matrix for
% the N-S equation that's a function of the 'u' velocity

F1C = zeros(N,M);

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        %% Top
        F1C(N,j) = -((u(N,j+1)*u(N,j+1) - u(N,j-1)*u(N,j-1))/(2*mesh.dx)) - ((2*u(N-1,j)*v(N-1,j))/(mesh.dy));
         
        % Left
        F1C(i,1) = -((u(i,2)*u(i,2)) - (-u(i,M-1)*-u(i,M-1)))/(2*mesh.dx) -((u(i+1,j)*v(i+1,j)) - (-u(i-1,j)*-v(i-1,j)))/(2*mesh.dy);             
        %% Main Domain
   
        F1C(i,j) = -((u(i,j+1)*u(i,j+1) - u(i,j-1)*u(i,j-1))/(2*mesh.dx)) - ((u(i+1,j)*v(i+1,j) - u(i-1,j)*v(i-1,j))/(2*mesh.dy));
        
        %% Bottom
        F1C(1,j) = 0;
        
        % Right
        F1C(i,M) = -((u(i,2)*u(i,2) - u(i,M-1)*u(i,M-1))/(2*mesh.dx)) - ((u(i+1,M)*v(i+1,M) - u(i-1,M)*v(i-1,M))/(2*mesh.dy));       
    end
end

% Top Left Corner
F1C(N,1) = -((u(N,2)*u(N,2) - u(N,M-1)*u(N,M-1))/(2*mesh.dx)) - ((2*u(N-1,1)*v(N-1,1))/(mesh.dy));

% Bottom Left Corner
F1C(1,1) = 0;

% Bottom Right Corner
F1C(1,M) = 0;

% Top Right Corner
F1C(N,M) = -((u(N,2)*u(N,2) - u(N,M-1)*u(N,M-1))/(2*mesh.dx)) - ((2*u(N-1,M))/(2*mesh.dy));
end       