function [F1V2]=F1V2(N,M,u,mesh,nonDimParams)

F1V2=zeros(M,N);

for r=1:M

        % bottom
        F1V2(r,1)=(1/nonDimParams.reynolds_)*((u(r,2)-u(r,2))/(mesh.dy)^2);
        
        % top
        F1V2(r,N)=(1/nonDimParams.reynolds_)*((u(r,N-1)-u(r,N)+u(r,N-1))/(mesh.dy)^2);
end 

for i=2:M-1
    for j=2:N-1
        
        % center
        F1V2(i,j)=(1/nonDimParams.reynolds_)*((u(i,j+1)-2*u(i,j)+u(i,j-1))/((mesh.dy)^2));   
    end 
end 
end 
