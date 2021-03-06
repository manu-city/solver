function [Pressure] = poisson(mesh, dt, pred, AP, b)
% Split each boundary into a seperate for loop for speed, change in
% orientation means that some of the '+M' become '-M', edited some of the b
% functions to incldue the boundary conditions instead of simply keeping
% them as adding or subtacting from eachother

M = mesh.nx;                         
N = mesh.ny;   
dx = mesh.dx;
dy = mesh.dy;
u = pred.VOFu;
v = pred.VOFv;

for i = 2:N-1     % ROWS
    % Centre of the domain
    for j = 2:M-1 % COLUMNS
        pointer = (i - 1)*M + j; 
        % b vector
        b(pointer,1) = 2/dt * ((u(i,j+1) - u(i,j-1))/(2*dx) + ...
                               (v(i+1,j) - v(i-1,j))/(2*dy));
    end
end

for j = M
    for i = 2:N-1
        % Right of domain
        pointer = (i - 1)*M + j;
        b(pointer,1) = 2/dt * ((u(i,2) - u(i,j-1))/(2*dx) + ...
                           (v(i+1,j) - v(i-1,j))/(2*dy));
    end
end

for j = 1
    for i = 2:N-1
        % Left of domain
        pointer = (i - 1)*M + j;
        b(pointer,1) = 2/dt * ((u(i,2) - u(i,M-1))/(2*dx) + ...
                           (v(i+1,j) - v(i-1,j))/(2*dy)); 
    end
end


for i = N
    for j = 2:M-1
        %Top of the domain
        pointer = (i - 1)*M + j;
        b(pointer,1) = 2/dt * ((u(i,j+1) - u(i,j-1))/(2*dx) + ...
                           (-v(N-1,j))/(dy));
    end
end

for i = 1
    for j = 2:M-1
        %Bottom of domain
        pointer = (i - 1)*M + j;
        b(pointer,1) = 2/dt * ((u(i,j+1) - u(i,j-1))/(2*dx) + ...
                           (+v(2,j))/(dy)); 
    end
end


for i = 1
    for j = M
        % Bottom Right Corner
        pointer = (i - 1)*M + j;
        % b vector
        b(pointer,1) = 0;
    end
end

for i = 1
    for j = 1
        % Bottom left corner
        pointer = (i - 1)*M + j;
        b(pointer,1) = 2/dt * ((u(1,2) - u(1,M-1))/(2*dx) + ...
                       (+v(2,1))/(dy));
    end
end

for i = N
    for j = M
        % Top right corner 
        pointer = (i - 1)*M + j;
        b(pointer,1) = 2/dt * ((u(N,2) - u(N,M-1))/(2*dx) + ...
                       (-v(N-1,M))/(dy));
    end
end

for i = N
    for j = 1
        % Top left corner
        pointer = (i - 1)*M + j;
        b(pointer,1) = 2/dt * ((u(1,2) - u(1,M-1))/(2*dx) + ...
                       (-v(N-1,1))/(dy));
    end
end

% Reform pressure matrix
Press = AP\b;
Pressure  = reshape(Press, [M,N])';

end