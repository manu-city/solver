function [pred] = obstacles_(mesh, nonDimParams, dimParams, pred)
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
dx = mesh.dx;
dy = mesh.dy;

o_x = nonDimParams.o_x_;
o_y = nonDimParams.o_y_;

o_num = dimParams.o_num;

u = pred.predU;
v = pred.predV;

if o_num == 0
    VOF = ones(N,M);
else
    % Total space between ribs
    spacing = floor(M - o_num * o_x);
    
    % Middle spaces
    mid_spacing = spacing/o_num;
    
    % height of the ribs
    o_height = floor(o_y/dy);
    
    % Width of ribs
    o_width = ceil(o_x/dx);
    
    % Obstacle space
    obs = zeros(length(1:1 + (o_height - 1)), length(floor(mid_spacing/2):ceil(M-(mid_spacing/2))));
    
    % space between ribs inside obstacle space
    space = ceil((size(obs,2) - o_width*o_num)/(o_num-1));
    
    % single space in between ribs
    dom_space = ones(size(obs,1), space + 1);
    
    for x = 1:o_num-1
        obs(:,o_width*x + (x-1)*space:x*(o_width + space)) = dom_space;
    end
    
    % domain
    VOF = ones(N,M);
    VOF(1:1 + (o_height - 1), floor(mid_spacing/2):ceil(M-(mid_spacing/2))) = obs;
end
pred.VOFu = u .* VOF;
pred.VOFv = v .* VOF;
end