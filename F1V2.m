function [F1V2]=F1V2(N,M,u,mesh,nonDimParams)

F1V2=zeros(N,M);
Re = nonDimParams.reynolds_;

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Main
        F1V2(i,j) = (1/Re)*((u(i+1,j) - 2*u(i,j) + u(i-1,j))/((mesh.dy)^2));
        
        % Left
        F1V2(i,1) = (1/Re)*((u(i+1,1) - 2*u(i,1) + u(i-1,1))/((mesh.dy)^2));
        
        % Right
        F1V2(i,M) = (1/Re)*((u(i+1,M) - 2*u(i,M) + u(i-1,M))/((mesh.dy)^2));
        
        % Bottom
        F1V2(N,j) = 0;
        
        % Top
        F1V2(1,j) = (1/Re)*((2*u(2,j) - 2*u(1,j))/((mesh.dy)^2));
    end 
end 
        % Bottom Left Corner
        F1V2(N,1) = 0;
        
        % Bottom Right Corner
        F1V2(N,M) = 0;
        
        % Top Left Corner
        F1V2(1,1) = (1/Re)*((2*u(2,1) - 2*u(1,1))/((mesh.dy)^2));
        
        % Top Right Corner
        F1V2(1,M) = (1/Re)*((2*u(2,M) - 2*u(1,M))/((mesh.dy)^2));
end 