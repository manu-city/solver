function [F1V1]=F1V1(N,M,u,mesh,nonDimParams)

F1V1=zeros(N,M);


for r=1:N
    
     
       % left
       F1V1(r,1)=(1/nonDimParams.reynolds_)*((u(r,2)-2*u(r,1)+u(r,M-1))/(mesh.dx)^2);
       
       % right
       F1V1(r,M)=(1/nonDimParams.reynolds_)*((u(r,2)-2*u(r,M)+u(r,M-1))/((mesh.dx)^2));
       
end


for i=2:N-1
    for j=2:M-1
        
        % center
        F1V1(i,j)=(1/nonDimParams.reynolds_)*((u(i,j+1)-2*u(i,j)+u(i,j-1))/((mesh.dx)^2));
     
    end 
end 
end 

