function [VOF, pred] = obstacles(mesh, nonDimParams, dimParams, N, M, pred)

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

for i = N-(o_height):N
    
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

pred.predU = u .* VOF;
pred.predV = v .* VOF;

end