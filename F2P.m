function [F2P] = F2P(N,M,P,mesh)
% Function file which creates the Divergence of Pressure Flux matrix for
% the N-S equation that's a function of the 'v' velocity

F2P = zeros(N,M);

for k = 1:M
    % Top
    F2P(1,k) = (-P(N-1,k) - P(N-1,k))/(2*mesh.dx);
    
    % Bottom
    F2P(N,k) = (-P(2,k) - P(2,k))/(2*mesh.dx);
end

for i = 2:N-1
    for j = 2:M-1
        % Main
        F2P(i,j) = (-P(i+1,j) - P(i-1,j))/(2*mesh.dy);
    end
end
