function [F2V2]=F2V2(N,M,v,mesh,nonDimParams)

F2V2=zeros(M,N);

for r=1:M

        %top
        F2V2(r,N)=(1/nonDimParams.reynolds_)*((v(r,N-1)+v(r,N-1))/(mesh.dy)^2);
        
        %bottom
        F2V2(r,1)=(1/nonDimParams.reynolds_)*((v(r,2)-v(r,2))/(mesh.dy)^2);
        
end 

for i=2:M-1
    for j=2:N-1
       
        % center
        F2V2(i,j)=(1/nonDimParams.reynolds_)*((v(i,j+1)-2*v(i,j)+v(i-1,j-1))/((mesh.dy)^2));
      
    end 
end 
