function [F1P] = F1P(N,M,P,mesh)
% Function file which creates the Divergence of Pressure Flux matrix for
% the N-S equation that's a function of the 'u' velocity

F1P = zeros(N,M);

for k = 1:N
    % Left
    F1P(k,1) = (-P(k,2) - P(k,N-1))/(2*mesh.dx);
    
    % Right
    F1P(k,M) = (-P(k,2) - P(k,N-1))/(2*mesh.dx);
end

for i = 2:N-1
    for j = 2:M-1
        P(i,j) = (-P(i+1,j) - P(i-1,j))/(2*mesh.dx);
    end
end