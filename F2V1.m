function [F2V1]=F2V1(N,M,v,mesh,nonDimParams)

F2V1=zeros(N,M);

for i=2:N-1
    for j=2:M-1
       
        %center
        F2V1(i,j)=(1/nonDimParams.reynolds_)*((v(i,j+1)-2*v(i,j)+v(i,j-1))/((mesh.dx)^2));
        
        %left
        F2V1(i,1)=(1/nonDimParams.reynolds_)*((v(i,2)-2*v(i,1)+v(i,M-1))/((mesh.dx)^2));
        
        %right
        F2V1(i,M)=(1/nonDimParams.reynolds_)*((v(i,2)-2*v(i,M)+v(i,M-1))/((mesh.dx)^2));
        
        
    end 
end 
