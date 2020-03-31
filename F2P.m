function [F2P] = F2P(N,M,P,mesh)
% Function file which creates the Divergence of Pressure Flux matrix for
% the N-S equation that's a function of the 'v' velocity

F2P = zeros(M,N);

for k = 1:M
    % Top
    F2P(k,N) = (-P(k,N-1) - P(k,N-1))/(2*mesh.dx);
    
    % Bottom
    F2P(k,1) = (-P(k,2) - P(k,2))/(2*mesh.dx);
end

for i = 2:M-1
    for j = 2:N-1
        % Main
        F2P(i,j) = (-P(i+1,j) - P(i-1,j))/(2*mesh.dy);
    end
end
