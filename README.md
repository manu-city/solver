# Prediction

## Dependencies
- mesh
- nonDimParams
- G
- W
- dt
- t
- k
- solution
- N
- M

## Initialisation of variables

    % Extracting Variables
    Re  = nonDimParams.reynolds_;
    dy  = mesh.dy;

    % Once first time step has passed
    if k > 1
        dt2 = t(k) - t(k-1);

        % Previous GP for time step > 1
        GP = G;
    else
        % Previous GP at time = 0
        GP = 0;
    end

    % Previous solution
    u_0 = solution.u;
    v_0 = solution.v;

    % Initialising predictions and A
    predU   = zeros(N,M);
    predV   = zeros(N,M);
    A       = zeros(N,1); % single coloumn vector due to A matrix generation 
                          % method below
                          
## Initial time step and A matrix beta coefficient 

    % initial time step
    if k == 1
        b1 = u_0 + 1.5 * dt * G.G1 + dt * W.W1;
        b2 = v_0 + 1.5 * dt * G.G2 + dt * W.W2;
    else
        b1 = u_0 + 1.5 * dt * G.G1 - 0.5 * dt2 * GP.G1 + dt * W.W1;
        b2 = v_0 + 1.5 * dt * G.G2 - 0.5 * dt2 * GP.G2 + dt * W.W2;
    end

    % beta coefficient
    beta = dt/(2*Re*dy^2);
    
## A matrix for u and v velocities

    % u-velocity A matrix
    A_u = full(spdiags([A - beta...
                        A + 1 + 2*beta...
                        A - beta], -1:1, N, N));

    % v-velocity A matrix
    A_v = A_u;

    % Top condition for u-velocity in A matrix 
    A_u(1, 1)         = (1 + 2*beta);
    A_u(1, 2)         = -2*beta;
    
## Prediction matrices calculation with boundary conditions

    for j = 1:M

        % Corrected x-velocity
        predU(:, j)   = A_u\b1(:, j);

        % Corrected y-velocity
        predV(:, j)   = A_v\b2(:, j);

        % Bottom (both u- and v-velocities = 0)
        predU(N, j) = 0;
        predV(N, j) = 0;

        % Top (v-velocity = 0)
        predV(1, j) = 0;

    end
    
## Storing

    prediction.predU = predU;
    prediction.predV = predV;


# Projection

## Dependencies
- mesh
- Pressure
- dt
- prediction
- N
- M

## Initialisation

        dx = mesh.dx;
        dy = mesh.dy;

        p = Pressure;
        pred_u = prediction.predU;
        pred_v = prediction.predV;

        u  = zeros(N,M);
        v  = zeros(N,M);

## Update of velocity in Main, right and left domain

        for i = 2:N-1

            % Main Domain
            for j = 2:M-1

                u(i,j) = pred_u(i,j) - dt/2 * (p(i,j+1) - p(i,j-1))/2*dx;
                v(i,j) = pred_v(i,j) - dt/2 * (p(i+1,j) - p(i-1,j))/2*dy;

            end

            % Right of domain
            j = M;
            u(i,j) = pred_u(i,j) - dt/2 * (p(i,2) - p(i,j-1))/2*dx;
            v(i,j) = pred_v(i,j) - dt/2 * (p(i+1,j) - p(i-1,j))/2*dy;

            % Left of domain
            j = 1;
            u(i,j) = pred_u(i,j) - dt/2 * (p(i,j+1) - p(i,M-1))/2*dx;
            v(i,j) = pred_v(i,j) - dt/2 * (p(i+1,j) - p(i-1,j))/2*dy;  

        end
        
 ## Update of velocity at the top and bottom domain, including corners
 
         for j = 2:M-1

            % Top of domain
            i = 1;
            u(i,j) = pred_u(i,j) - dt/2 * (p(i,j+1) - p(i,j-1))/2*dx;
            v(i,j) = 0;

            % Bottom of domain
            i = N;
            u(i,j) = 0;
            v(i,j) = 0;

        end

        % Top Right
        i = 1;
        j = M;
        u(i,j) = pred_u(i,j) - dt/2 * (p(i,2) - p(i,j-1))/2*dx;

        % Top Left
        j = 1;
        u(i,j) = pred_u(i,j) - dt/2 * (p(i,j+1) - p(i,M-1))/2*dx;
        
## Storing solution

        solution.u = u;
        solution.v = v;
        solution.P = p;
 
 
