function [prediction] = prediction(mesh, nonDimParams, G, W, dt, t, k, solution, N, M)

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
                      
% b matrix

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

% u-velocity A matrix
A_u = full(spdiags([A - beta ...
                    A + 1 + 2*beta ...
                    A - beta], -1:1, N, N));
           
% v-velocity A matrix
A_v = A_u;

% Top condition for u-velocity in A matrix 
A_u(1, 1)         = (1 + 2*beta);
A_u(1, 2)         = -2*beta;

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

prediction.predU = predU;
prediction.predV = predV;
end