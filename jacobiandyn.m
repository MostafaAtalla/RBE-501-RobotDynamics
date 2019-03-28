function J=jacobiandyn(Tdyn,CoM,index)
%First get the origin of the tip frame
DOF=size(Tdyn,2);
J=sym(zeros(3,DOF));
Ocom=sym(zeros(3,DOF));
o=sym(zeros(3,DOF));
z=sym(zeros(3,DOF));
z(:,1)=sym([0;0;1]);
for i=1:DOF
    k=(Tdyn{index,i}*[CoM(1,i);0;0;1]);
    Ocom(:,i)=k(1:3,1);
end
for i=1:(DOF-1)
z(:,i+1)=Tdyn{index,i}(1:3,3);
o(:,i+1)=Tdyn{index,i}(1:3,4);
end
%% Computing the Jacobian Matrix %%
for i=1:DOF
    x=z(:,i);
    y=simplify(Ocom(:,i)-o(:,i));
    J(1:3,i)=simplify(cross(x,y));
end
end

