function [F1C] = F1C(N,M,u,v,mesh)
% Function file which creates the Divergence of Convective Flux matrix for
% the N-S equation that's a function of the 'u' velocity

F1C = zeros(M,N);

for i = 2:M-1  
    for j = 2:N-1
        %% Top
        F1C(i,N) = -((u(i+1,N)*u(i+1,N) - u(i-1,N)*u(i-1,N))/(2*mesh.dx)) - ((u(i,N-1)*v(i,N-1))/2*mesh.dy); 
        
        % Top Left Corner
        F1C(1,N) = -((u(2,N)*u(2,N) - u(M-1,N)*u(M-1,N))/(2*mesh.dx)) - ((u(1,N-1)*v(1,N-1))/(2*mesh.dy));
        
        % Left
        F1C(1,j) = -((u(2,j)*u(2,j) - u(M-1,j)*u(M-1,j))/(2*mesh.dx)) - ((u(1,j+1)*v(1,j+1))/(2*mesh.dy));
        
        % Bottom Left Corner
        F1C(1,1) = -((u(1,2)*v(1,2) - (-u(1,2)*-v(1,2)))/(2*mesh.dy));
        
        %% Main Domain
        F1C(i,j) = -((u(i+1,j)*u(i+1,j) - u(i-1,j)*u(i-1,j))/(2*mesh.dx)) - ((u(i,j+1)*v(i,j+1) - u(i,j-1)*v(i,j-1))/(2*mesh.dy));
        
        %% Bottom
        F1C(i,1) = -((u(i,2)*v(i,2)) - (-u(i,2)*-v(i,2)))/(2*mesh.dy);
        
        % Bottom Right Corner
        F1C(M,1) = -((u(M,2)*v(M,2)) - (-u(M,2)*-v(M,2)))/(2*mesh.dy);
        
        % Right
        F1C(M,j) = -((u(2,j)*u(2,j)) - (u(M-1,j)*u(M-1,j))/(2*mesh.dx)) - ((u(M,j+1)*v(M,j+1)) - (-u(M,j-1)*v(M,j-1)))/(2*mesh.dy);
         
        % Top Right Corner
        F1C(M,N) = -((-u(2,N)*u(2,N)) - (u(M-1,N)*u(M-1,N))/(2*mesh.dx)) - ((u(M,N-1)*v(M,N-1))/(2*mesh.dy));
        
    end
end