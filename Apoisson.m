function [AP] = Apoisson(mesh, A)

M = mesh.nx;                         
N = mesh.ny;   
dx = mesh.dx;
dy = mesh.dy;

for i = 2:N-1     % ROWS
    % Centre of the domain
    for j = 2:M-1 % COLUMNS
        pointer = (i - 1)*M + j;
        % A matrix
        A(pointer, pointer - M) = 1/dy^2;
        A(pointer, pointer - 1) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer + 1) = 1/dx^2;
        A(pointer, pointer + M) = 1/dy^2;  
    end
end

for j = M
    for i = 2:N-1
        % Right of domain
        pointer = (i - 1)*M + j;
        A(pointer, pointer - M) = 1/dy^2;
        A(pointer, pointer - (M-2)) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer - 1) = 1/dx^2;
        A(pointer, pointer + M) = 1/dy^2;
    end
end

for j = 1
    for i = 2:N-1
        % Left of domain
        pointer = (i - 1)*M + j;
        A(pointer, pointer - M) = 1/dy^2;
        A(pointer, pointer + 1) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer + (M-2)) = 1/dx^2;  
        A(pointer, pointer + M) = 1/dy^2;
    end
end


for i = N
    for j = 2:M-1
        %Top of the domain
        pointer = (i - 1)*M + j;
        A(pointer, pointer - M) = 2/dy^2;
        A(pointer, pointer - 1) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer + 1) = 1/dx^2;
    end
end

for i = 1
    for j = 2:M-1
        %Bottom of domain
        pointer = (i - 1)*M + j;
        A(pointer, pointer + M) = 2/dy^2;
        A(pointer, pointer - 1) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer + 1) = 1/dx^2;  
    end
end


for i = 1
    for j = M
        % Bottom Right Corner
        pointer = (i - 1)*M + j;
        % A matrix
        A(pointer, pointer) = 1;
    end
end

for i = 1
    for j = 1
        % Bottom left corner
        pointer = (i - 1)*M + j;
        A(pointer, pointer + M) = 2/dy^2;
        A(pointer, pointer + 1) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer + (M-2)) = 1/dx^2;
    end
end

for i = N
    for j = M
        % Top right corner 
        pointer = (i - 1)*M + j;
        A(pointer, pointer - M) = 2/dy^2;
        A(pointer, pointer - (M-2)) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer - 1) = 1/dx^2;
    end
end

for i = N
    for j = 1
        % Top left corner
        pointer = (i - 1)*M + j;
        A(pointer, pointer - M) = 2/dy^2;
        A(pointer,pointer + 1) = 1/dx^2;
        A(pointer,pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer,pointer + (M-2)) = 1/dx^2;
    end
end

% Sparse the matrix for speed
AP = sparse(A);
end

