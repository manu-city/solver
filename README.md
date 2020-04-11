# Poisson Solver

## Dependencies
- mesh
- N
- M
- dt
- solution

## Initialisation

    dx = mesh.dx;
    dy = mesh.dy;
    A = zeros(N*M);
    b = zeros(N*M, 1);
    u = solution.u;
    v = solution.v;
    
## Centre, Right and Left domain A matrix and b vector

    for i = 2:N-1
        % Centre of the domain
        for j = 2:M-1
            pointer = (i - 1)*M + j;
            % A matrix
            A(pointer, pointer - M) = 1/dy^2;
            A(pointer, pointer - 1) = 1/dx^2;
            A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
            A(pointer, pointer + 1) = 1/dx^2;
            A(pointer, pointer + M) = 1/dy^2;  
            % b vector
            b(pointer,1) = 2/dt * ((u(i,j+1) - u(i,j-1))/2*dx + ...
                                   (v(i+1,j) - v(i-1,j))/2*dy);
        end

        % Right of domain
        j = M;
        pointer = (i - 1)*M + j;
        A(pointer, pointer - M) = 1/dy^2;
        A(pointer, pointer - (M-2)) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer - 1) = 1/dx^2;
        A(pointer, pointer + M) = 1/dy^2;
        b(pointer,1) = 2/dt * ((u(i,2) - u(i,j-1))/2*dx + ...
                               (v(i+1,M) - v(i-1,M))/2*dy);

        % Left of domain
        j = 1;
        pointer = (i - 1)*M + j;
        A(pointer, pointer - M) = 2/dy^2;
        A(pointer, pointer + 1) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer + (M-2)) = 1/dx^2;  
        A(pointer, pointer + M) = 1/dy^2;
        b(pointer,1) = 2/dt * ((u(i,2) - u(i,M-1))/2*dx + ...
                               (v(i+1,1) - v(i-1,1))/2*dy); 
    end
    
## Top and bottom domain A matrix and b vector

    for j = 2:M-1
        %Top of the domain
        i = 1;
        pointer = (i - 1)*M + j;
        A(pointer, pointer + M) = 2/dy^2;
        A(pointer, pointer - 1) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer + 1) = 1/dx^2;
        b(pointer,1) = 2/dt * ((u(i,j+1) - u(i,j-1))/2*dx + ...
                               (v(2,j) + v(2,j))/2*dy);

        %Bottom of domain
        i = N;
        pointer = (i - 1)*M + j;
        A(pointer, pointer - M) = 2/dy^2;
        A(pointer, pointer - 1) = 1/dx^2;
        A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
        A(pointer, pointer + 1) = 1/dx^2;  
        b(pointer,1) = 2/dt * ((u(N,j+1) - u(N,j-1))/2*dx + ...
                               (-v(N-1,j) - v(N-1,j))/2*dy); 
    end

## All 4 domain corners A matrix and b vector

    for i = [N, 1] 
        for j = [M, 1]
            ij = i+j;
            switch ij 
                % Bottom right corner
                case N + M
                    pointer = (i - 1)*M + j;
                    % A matrix
                    A(pointer, pointer - M) = 1/dy^2;
                    A(pointer, pointer - (M-2)) = 1/dx^2;
                    A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
                    A(pointer, pointer - 1) = 1/dx^2;
                    % b vector
                    b(pointer,1) = 2/dt * ((u(N,2) - u(N,M-1))/2*dx + ...
                                           (-v(N-1,M) - v(N-1,M))/2*dy);

                % Bottom left corner
                case N + 1
                    pointer = (i - 1)*M + j;
                    A(pointer, pointer - M) = 2/dy^2;
                    A(pointer, pointer + 1) = 1/dx^2;
                    A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
                    A(pointer, pointer + (M-2)) = 1/dx^2;
                    b(pointer,1) = 2/dt * ((u(N,2) - u(N,M-1))/2*dx + ...
                                           (v(N-1,1) - v(N-1,1))/2*dy);

                % Top right corner 
                case 1 + M
                    pointer = (i - 1)*M + j;
                    A(pointer, pointer + M) = 2/dy^2;
                    A(pointer, pointer - (M-2)) = 1/dx^2;
                    A(pointer, pointer) = -(2/dx^2 + 2/dy^2);
                    A(pointer, pointer - 1) = 1/dx^2;
                    b(pointer,1) = 2/dt * ((u(1,2) - u(1,M-1))/2*dx + ...
                                           (v(2,M) + v(2,M))/2*dy);

                % Top left corner   
                case 1 + 1
                    pointer = (i - 1)*M + j;
                    A(pointer, pointer) = 1;
                    b(pointer,1) = 0
            end
        end
    end
    
## Sparse A
    A = sparse(A);
    
## Reform pressure matrix
    Pressure = reshape(A\b, [M,N])';
  
