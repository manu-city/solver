function [pred] = obstacles(mesh, nonDimParams, dimParams, pred)
M = mesh.nx;   % Gridpoints in x-axis (Columns)                       
N = mesh.ny;   % Gridpoints in y-axis (Rows)
dx = mesh.dx;
dy = mesh.dy;

o_x = nonDimParams.o_x_;
o_y = nonDimParams.o_y_;
clearance = nonDimParams.clear_;

o_num = dimParams.o_num;

u = pred.predU;
v = pred.predV;

% Grid points in x in between ribs
spacing = clearance/dx;

% Grid points in x rib occupies
o_length = o_x/dx;

% Grid points in y rib occupies
o_height = o_y/dy;

% domain
VOF = ones(N,M);

for i = 1:1+(o_height)
    for n = 1:o_num
        switch n
            case 1 
                for j = spacing/2 + 1:(spacing/2 + 1) + o_length
                    VOF(i,j) = 0;    
                end
            otherwise 
                for j = (spacing + o_length) * (n) - (spacing/2 + o_length - 1) : (spacing + o_length) * (n) - (spacing/2 + o_length - 1) + o_length
                    VOF(i,j) = 0;
                end
        end 
    end
end

pred.VOFu = u .* VOF;
pred.VOFv = v .* VOF;
end