function [F1V1]=F1V1(N,M,u,mesh,nonDimParams)

F1V1=zeros(M,N);


for r=2:N
    
     
      % left
       F1V1(1,r)=(1/nonDimParams.reynolds_)*((u(2,r)-2*u(1,r)+u(N-1,r))/((mesh.dx)^2));
       
       % right
       F1V1(M,r)=(1/nonDimParams.reynolds_)*((u(2,r)-2*u(N,r)+u(N-1,r))/((mesh.dx)^2));
       
end


for i=2:M-1
    for j=2:N-1
        
        % center
        F1V1(i,j)=(1/nonDimParams.reynolds_)*((u(i+1,j)-2*u(i,j)+u(i-1,j))/((mesh.dx)^2));
    
    end 
end 
end 

