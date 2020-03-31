function [F1P] = F1P(N,M,P,mesh)
% Function file which creates the Divergence of Pressure Flux matrix for
% the N-S equation that's a function of the 'u' velocity

F1P = zeros(M,N);

for k = 1:N
    % Left
    F1P(1,k) = (-P(2,k) - P(N-1,k))/(2*mesh.dx);
    
    % Right
    F1P(M,k) = (-P(2,k) - P(N-1,k))/(2*mesh.dx);
end

for i = 2:M-1
    for j = 2:N-1
        P(i,j) = (-P(i+1,j) - P(i-1,j))/(2*mesh.dx);
    end
end