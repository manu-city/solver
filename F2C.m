function [F2C] = F2C(N,M,u,v,mesh)
% Function file which creates the Divergence of Convective Flux matrix for
% the N-S equation that's a function of the 'v' velocity

F2C = zeros(N,M);

for i = 2:N-1
    for j = 2:M-1
        %% Top
        F2C(1,j) = - ((v(2,j)*v(2,j) - v(2,j)*v(2,j))/(2*mesh.dy));
        
        % Top Left Corner
        F2C(1,1) = -((v(2,1)*v(2,1) - (v(2,1)*v(2,1)))/(2*mesh.dy));
        
        % Left
        F2C(i,1) = -((u(i,2)*v(i,2) - (u(i,M-1)*v(i,M-1)))/(2*mesh.dx)) - (v(i+1,1)*v(i+1,1) - v(i-1,1)*v(i-1,1))/(2*mesh.dy);
        
        % Bottom Left Corner
        F2C(N,1) = -((v(N-1,1)*v(N-1,1)) - (v(N-1,1)*v(N-1,1)))/(2*mesh.dy);
        
        %% Main Domain
        F2C(i,j) = ((u(i,j+1)*v(i,j+1) - u(i,j+1)*v(i,j+1))/(2*mesh.dx)) - ((v(i+1,j)*v(i+1,j) - v(i-1,j)*v(i-1,j))/(2*mesh.dy));
        
        %% Bottom
        F2C(N,j) = - ((v(N-1,j)*v(N-1,j) - v(N-1,j)*v(N-1,j))/(2*mesh.dy));
        
        % Bottom Right Corner
        F2C(N,M) = - ((v(N-1,M)*v(N-1,M) - v(N-1,M)*v(N-1,M))/(2*mesh.dy));
        
        % Right
        F2C(i,M) = -(u(i,2)*v(i,2) - u(i,M-1)*v(i,M-1))/(2*mesh.dx) - (v(i+1,M)*v(i+1,M) - v(i-1,M)*v(i-1,M))/(2*mesh.dy);
        
        % Top Right Corner
        F2C(1,M) = - ((v(2,j)*v(2,j) - v(2,j)*v(2,j))/(2*mesh.dy));
        
    end
end