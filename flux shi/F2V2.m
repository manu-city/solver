function [F2V2]=F2V2(N,M,v,mesh,nonDimParams)

F2V2=zeros(N,M);
Re = nonDimParams.reynolds_;

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Center
        F2V2(i,j) = (1/Re)*((v(i+1,j)-2*v(i,j)+v(i-1,j))/((mesh.dy)^2));
        
        % Left
        F2V2(i,1) = (1/Re)*((v(i+1,1)-2*v(i,1)+v(i-1,1))/((mesh.dy)^2));
        
        % Right
        F2V2(i,M) = (1/Re)*((v(i+1,M)-2*v(i,M)+v(i-1,M))/((mesh.dy)^2));
        
        % Bottom
        F2V2(N,j) = (1/Re)*((v(N-1,j)-2*v(N,j))/((mesh.dy)^2));
        
        % Top
        F2V2(1,j) = (1/Re)*((2*v(2,j))/((mesh.dy)^2));
    end 
end 
        % Bottom Left Corner
        F2V2(N,1) = (1/Re)*((v(N-1,1)-2*v(N,1))/((mesh.dy)^2));
        
        % Bottom Right Corner
        F2V2(N,M) = (1/Re)*((v(N-1,M)-2*v(N,M))/((mesh.dy)^2));
        
        % Top Left Corner
        F2V2(1,1) = (1/Re)*((2*v(2,1))/((mesh.dy)^2));
        
        % Top Right Corner
        F2V2(1,M) = (1/Re)*((2*v(2,M))/((mesh.dy)^2));
end 