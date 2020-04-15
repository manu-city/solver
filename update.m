function [solution] = update(mesh, Pressure, dt, prediction, N, M)

dx = mesh.dx;
dy = mesh.dy;

p = Pressure;
pred_u = prediction.predU;
pred_v = prediction.predV;

u  = zeros(N,M);
v  = zeros(N,M);

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
    i = 1;
    u(i,j) = pred_u(i,j) - (dt/2) * (p(i,j+1) - p(i,j-1))/(2*dx);
    v(i,j) = 0;
    
    
end

% Top Right
i = 1;
j = M;
u(i,j) = pred_u(i,j) - (dt/2) * (p(i,2) - p(i,j-1))/(2*dx);

% Top Left
j = 1;
u(i,j) = pred_u(i,j) - (dt/2) * (p(i,j+1) - p(i,M-1))/(2*dx);

% Bottom of domain
i = N;
u(i,:) = 0;
v(i,:) = 0;

solution.u = u;
solution.v = v;
solution.P = p;

end