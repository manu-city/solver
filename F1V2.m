function [F1V2]=F1V2(N,M,u,mesh,nonDimParams)

F1V2=zeros(N,M);

for r=1:M
        % bottom
        F1V2(N,r)=(1/nonDimParams.reynolds_)*((u(N-1,r) + u(N-1,r))/(mesh.dy)^2);
        % top
        F1V2(1,r)=(1/nonDimParams.reynolds_)*((u(2,r)-u(1,r)+ u(2,r))/(mesh.dy)^2);
end 

for i=2:N-1
    for j=2:M-1     
        % center
        F1V2(i,j)=(1/nonDimParams.reynolds_)*((u(i+1,j)-2*u(i,j)+u(i-1,j))/((mesh.dy)^2));   
    end 
end 
end 
