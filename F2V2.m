function [F2V2]=F2V2(N,M,v,mesh,nonDimParams)

F2V2=zeros(N,M);

for r=1:M

        %top
        F2V2(1,r)=(1/nonDimParams.reynolds_)*((v(2,r)+v(2,r))/(mesh.dy)^2);
        
        %bottom
        F2V2(N,r)=(1/nonDimParams.reynolds_)*((v(N-1,r)-v(N-1,r))/(mesh.dy)^2);
        
end 

for i=2:N-1
    for j=2:M-1
       
        % center
        F2V2(i,j)=(1/nonDimParams.reynolds_)*((v(i+1,j)-2*v(i,j)+v(i-1,j))/((mesh.dy)^2));
      
    end 
end 
