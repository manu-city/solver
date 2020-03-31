function [F2C] = F2C(N,M,u,v,mesh)
% Function file which creates the Divergence of Convective Flux matrix for
% the N-S equation that's a function of the 'v' velocity

F2C = zeros(M,N);

for i = 2:M-1
    for j = 2:N-1
        %% Top
        F2C(i,N) = -(v(i,N-1)*v(i,N-1)/(2*mesh.dy));
        
        % Top Left Corner
        F2C(1,N) = -(v(1,N-1)*v(1,N-1)/(2*mesh.dy));
        
        % Left
        F2C(1,j) = -((u(2,j)*v(2,j) - (u(M-1,j)*v(M-1,j)))/(2*mesh.dx)) - ((v(1,j+1)*v(1,j+1) - v(1,j-1)*v(1,j-1))/(2*mesh.dy));
        
        % Bottom Left Corner
        F2C(1,1) = ((v(1,2)*v(1,2) - (-v(1,2)*-v(1,2)))/(2*mesh.dy));
        
        %% Main Domain
        F2C(i,j) = ((u(i+1,j)*v(i+1,j) - u(i-1,j)*v(i-1,j))/(2*mesh.dx)) - ((v(i,j+1)*v(i,j+1) - v(i,j-1)*v(i,j-1))/(2*mesh.dy));
        
        %% Bottom
        F2C(i,1) = -((u(i+1,1)*v(i+1,1) - (u(i-1,1)*v(i-1,1)))/(2*mesh.dx)) - ((v(i,2)*v(i,2))/(2*mesh.dy));
        
        % Bottom Right Corner
        F2C(M,1) = ((v(M,2)*v(M,2)) - (-v(M,2)*-v(M,2)))/(2*mesh.dy);
        
        % Right
        F2C(M,j) = -((u(2,j)*v(2,j)) - (u(M-1,j)*v(M-1,j)))/(2*mesh.dx) - ((v(M,j+1)*v(M,j+1) - v(M,j-1)*v(M,j-1))/(2*mesh.dy));
        
        % Top Right Corner
        F2C(M,N) = (-v(M,N-1)*v(M,N-1))/(2*mesh.dy);
        
    end
end