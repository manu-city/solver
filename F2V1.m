function [F2V1]=F2V1(N,M,v,mesh,nonDimParams)

F2V1=zeros(M,N);

for i=2:M-1
    for j=2:N-1
       
        %center
        F2V1(i,j)=(1/nonDimParams.reynolds)*((v(i+1,j)-2*v(i,j)+v(i-1,j))/((mesh.dx)^2));
        
        %left
        F2V1(1,j)=(1/nonDimParams.reynolds)*((v(2,j)-2*v(1,j)+v(N-1,j))/((mesh.dx)^2));
        
        %right
        F2V1(M,j)=(1/nonDimParams.reynolds)*((v(2,j)-2*v(N,j)+v(N-1,j))/((mesh.dx)^2));
        
        
    end 
end 
