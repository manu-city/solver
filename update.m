function [solution] = update(mesh, Pressure, dt, pred, solution)
% Switched the neccesary boundaries due to the orientation
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
dx = mesh.dx;
dy = mesh.dy;

p = Pressure;
pred_u = pred.VOFu; % Used the VOF*u_star velocity instead of just the u_star
pred_v = pred.VOFv; % Used the VOF*v_star velocity instead of just the v_star

u  = solution.u;
v  = solution.v;

for i = 2:N-1       % ROWS
    
    % Main Domain
    for j = 2:M-1   % COLUMNS
        
        u(i,j) = pred_u(i,j) - (dt/2) * (p(i,j+1) - p(i,j-1))/(2*dx);
        v(i,j) = pred_v(i,j) - (dt/2) * (p(i+1,j) - p(i-1,j))/(2*dy);
        
    end
    
    % Right of domain
    j = M;
    u(i,j) = pred_u(i,j) - (dt/2) * (p(i,2) - p(i,j-1))/(2*dx);
    v(i,j) = pred_v(i,j) - (dt/2) * (p(i+1,j) - p(i-1,j))/(2*dy);
    
    % Left of domain
    j = 1;
    u(i,j) = pred_u(i,j) - (dt/2) * (p(i,j+1) - p(i,M-1))/(2*dx);
    v(i,j) = pred_v(i,j) - (dt/2) * (p(i+1,j) - p(i-1,j))/(2*dy);  
    
end

for j = 2:M-1
    
    % Top of domain
    i = N;
    u(i,j) = pred_u(i,j) - (dt/2) * (p(i,j+1) - p(i,j-1))/(2*dx);
    v(i,j) = 0;
    
    
end

% Top Right
i = N;
j = M;
u(i,j) = pred_u(i,j) - (dt/2) * (p(i,2) - p(i,j-1))/(2*dx);

% Top Left
i = N;
j = 1;
u(i,j) = pred_u(i,j) - (dt/2) * (p(i,j+1) - p(i,M-1))/(2*dx);

% Bottom of domain
for i = 1
    for j = 1:M
        u(i,j) = 0;
        v(i,j) = 0;
    end
end
solution.u = u;
solution.v = v;
solution.P = p;
end