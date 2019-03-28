function J=jacobiancom(T,index)
%First get the origin of the tip frame
DOF=size(T,2);
J=sym(zeros(6,DOF));
o=sym(zeros(3,DOF));
z=sym(zeros(3,DOF));
z(:,1)=sym([0;0;1]);
otip=T{index,DOF}(1:3,4);
for i=1:(DOF-1)
z(:,i+1)=T{index,i}(1:3,3);
o(:,i+1)=T{index,i}(1:3,4);
end
%% Computing the Jacobian Matrix %%
for i=1:DOF
    x=z(:,i);
    y=simplify(otip-o(:,i));
    J(1:3,i)=simplify(cross(x,y));
    J(4:6,i)=x;
end
end

