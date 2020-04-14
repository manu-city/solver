function [F2P] = F2P(N,M,P,mesh)
% Function file which creates the Divergence of Pressure Flux matrix for
% the N-S equation that's a function of the 'u' velocity

F2P = zeros(N,M);

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Main
        F2P(i,j) = -(P(i+1,j) - P(i-1,j))/(2*mesh.dy);
        
        % Bottom
        F2P(1,j) = 0;
        
        % Top
        F2P(N,j) = 0;
        
        % Left
        F2P(i,1) = -(P(i,2) - P(i,M-1))/(2*mesh.dy);
        
        % Right
        F2P(i,M) = -(P(i,2) - P(i,M-1))/(2*mesh.dy);
    end
end
        % Bottom Left
        F2P(1,1) = 0;
        
        % Bottom Right
        F2P(1,M) = 0;
        
        % Top Left
        F2P(N,1) = 0;
        
        % Top Right
        F2P(N,M) = 0;
end