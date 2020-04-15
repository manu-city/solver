function [VOF, prediction] = obstacles(mesh, nonDimParams, dimParams, N, M, prediction)

dx = mesh.dx;
dy = mesh.dy;

o_x = nonDimParams.o_x_;
o_y = nonDimParams.o_y_;
clear = nonDimParams.clear_;

o_num = dimParams.o_num;

u = prediction.predU;
v = prediction.predV;

% Grid points in x in between ribs
spacing = clear/dx;

% Grid points in x rib occupies
o_length = o_x/dx;

% Grid points in y rib occupies
o_height = o_y/dy;

% domain
VOF = ones(N,M);

for i = N-(o_height - 1):N
    
    for n = 1:o_num
        switch n
            case 1 
                for j = spacing/2 + 1:(spacing/2 + 1) + o_length
                    VOF(i,j) = 0;    
                end
            otherwise 
                for j = spacing * (n) + (2*n-3):spacing * (n) + (2*n-3) + o_length
                    VOF(i,j) = 0;
                end
        end 
    end
end

prediction.predU = u .* VOF;
prediction.predV = v .* VOF;

end