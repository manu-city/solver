function [F1V1]=F1V1(N,M,u,mesh,nonDimParams)

F1V1=zeros(N,M);
Re = nonDimParams.reynolds_;

for j = 2:M-1         % COLUMNS
    for i = 2:N-1     % ROWS
        % Center
        F1V1(i,j) = (1/Re)*((u(i,j+1) - 2*u(i,j) + u(i,j-1))/((mesh.dx)^2));
        
        % Left
        F1V1(i,1) = (1/Re)*((u(i,2) - 2*u(i,1) + u(i,M-1))/((mesh.dx)^2));
        
        % Right
        F1V1(i,M) = (1/Re)*((u(i,2) - 2*u(i,M) + u(i,M-1))/((mesh.dx)^2));
        
        % Bottom
        F1V1(N,j) = 0;
        
        % Top
        F1V1(1,j) = (1/Re)*((u(1,j+1) - 2*u(1,j) + u(1,j-1))/((mesh.dx)^2));
    end 
end 
        % Bottom Left Corner
        F1V1(N,1) = 0;
        
        % Bottom Right Corner
        F1V1(N,M) = 0;
        
        % Top Left Corner
        F1V1(1,1) = (1/Re)*((u(1,2) - 2*u(1,1) + u(1,M-1))/((mesh.dx)^2));
        
        % Top Right Corner
        F1V1(1,M) = (1/Re)*((u(1,2) - 2*u(1,M) + u(1,M-1))/((mesh.dx)^2));
end 