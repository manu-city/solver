function [F1C] = F1C(N,M,u,v,mesh)
% Function file which creates the Divergence of Convective Flux matrix for
% the N-S equation that's a function of the 'u' velocity

F1C = zeros(N,M);

for i = 2:N-1  
    for j = 2:M-1
        %% Top
        F1C(1,j) = -((u(2,j)*u(2,j) - u(N-1,j)*u(N-1,j))/(2*mesh.dx)) - ((u(1,j+1)*v(1,j+1))/(2*mesh.dy));
         
        % Top Left Corner
        F1C(1,1) = -((u(1,2)*v(1,2) - (-u(1,2)*-v(1,2)))/(2*mesh.dy));
         
        % Left
        F1C(i,1) = -((u(i,2)*v(i,2)) - (-u(i,2)*-v(i,2)))/(2*mesh.dy);
        
        % Bottom Left Corner
        F1C(N,1) = -((u(N,2)*v(N,2)) - (-u(N,2)*-v(N,2)))/(2*mesh.dy);
        
        %% Main Domain
        F1C(i,j) = -((u(i+1,j)*u(i+1,j) - u(i-1,j)*u(i-1,j))/(2*mesh.dx)) - ((u(i,j+1)*v(i,j+1) - u(i,j-1)*v(i,j-1))/(2*mesh.dy));
        
        %% Bottom
        F1C(N,j) = -((u(2,j)*u(2,j)) - (u(N-1,j)*u(N-1,j))/(2*mesh.dx)) - ((u(N,j+1)*v(N,j+1)) - (-u(N,j-1)*v(N,j-1)))/(2*mesh.dy);
        
        % Bottom Right Corner
        F1C(N,M) = -((-u(2,M)*u(2,M)) - (u(N-1,M)*u(N-1,M))/(2*mesh.dx)) - ((u(N,M-1)*v(N,M-1))/(2*mesh.dy));
        
        % Right
        F1C(i,M) = -((u(i+1,M)*u(i+1,M) - u(i-1,M)*u(i-1,M))/(2*mesh.dx)) - ((u(i,M-1)*v(i,M-1))/2*mesh.dy); 
        
        % Top Right Corner
        F1C(1,M) = -((u(2,M)*u(2,M) - u(N-1,M)*u(N-1,M))/(2*mesh.dx)) - ((u(1,M-1)*v(1,M-1))/(2*mesh.dy));
        
    end
end