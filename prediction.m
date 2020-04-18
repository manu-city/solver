function [pred] = prediction(mesh, nonDimParams, G, W, dt, t, k, solution, pred, predMat)
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
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
predU = pred.predU;
predV = pred.predV;
A1 = predMat.A1;
A1_ = predMat.A1_;
                    
% b matrix

% initial time step
if k == 1
    b1 = u_0 + dt * G.G1 + dt * W.W1;
    b2 = v_0 + dt * G.G2 + dt * W.W2;
else
    b1 = u_0 + 1.5 * dt * G.G1 - 0.5 * dt2 * GP.G1 + dt * W.W1;
    b2 = v_0 + 1.5 * dt * G.G2 - 0.5 * dt2 * GP.G2 + dt * W.W2;
end

% beta coefficient
beta = dt/(2*Re*dy^2);

% u-velocity A matrix
A_u = full(spdiags([A1 - beta ...
                    A1 + 1 + 2*beta ...
                    A1 - beta], -1:1, N-1, N-1));
           
% v-velocity A matrix
A_v = full(spdiags([A1_ - beta ...
                    A1_ + 1 + 2*beta ...
                    A1_ - beta], -1:1, N-2, N-2));

% Top condition for u-velocity in A matrix 
% Since the orienetation has changed the neumann boundary condition of the
% coefficient would be at the bottom of the A-Matrix now
A_u(N-1,N-1)         = (1 + 2*beta);
A_u(N-1,N-2)         = -2*beta;


for j = 1:M
    
    % Corrected x-velocity
    predU(2:N, j)   = A_u\b1(2:N, j);
    
    % Corrected y-velocity
    predV(2:N-1, j)   = A_v\b2(2:N-1, j);
end

pred.predU = predU;
pred.predV = predV;
end