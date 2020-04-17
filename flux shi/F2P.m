function [F2P] = F2P(N,M,P,mesh)
% Function file which creates the Divergence of Pressure Flux matrix for
% the N-S equation that's a function of the 'u' velocity

F2P = zeros(N,M);

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Main
        F2P(i,j) = -(P(i+1,j) - P(i-1,j))/(2*mesh.dy);
        
        % Bottom
        F2P(N,j) = -(- P(N-1,j))/(mesh.dy);
        
        % Top
        F2P(1,j) = -(P(2,j))/(mesh.dy);
        
        % Left
        F2P(i,1) = -(P(i+1,1) - P(i-1,1))/(2*mesh.dy);
        
        % Right
        F2P(i,M) = -(P(i+1,M) - P(i-1,M))/(2*mesh.dy);
    end
end
        % Bottom Left
        F2P(N,1) = -(- P(N-1,1))/(mesh.dy);
        
        % Bottom Right
        F2P(N,M) = -(- P(N-1,M))/(mesh.dy);
        
        % Top Left
        F2P(1,1) = -(P(2,1))/(mesh.dy);
        
        % Top Right
        F2P(1,M) = -(P(2,M))/(mesh.dy);
end