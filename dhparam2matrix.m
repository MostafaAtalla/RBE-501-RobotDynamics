% Function to compute the Transformation matrices based on the DH
% parameters %
% The inputs are the DH paramters theta,d,a,alpha
% The outputs are the transformation between each frame and the composite
% transfromation matrix.
% This function is made for HW2 RBE501 
% made by Mostafa Atalla
function [T] = dhparam2matrix(theta,d,a,alpha)
DOF=size(theta,2); %Compute the DOF based on the DH vectors size

%% Computing the Homogenous Transformation of the robot %%
T=cell(1,DOF);      %Defining homogeneous transformation cell to contation HTs

for i=1:DOF
ca=cos(alpha(1,i)); %Defining the cosine of the link twist
sa=sin(alpha(1,i)); %Defining the sine of the link twist
ct=cos(theta(1,i)); %Defining the cosine of the joint angle
st=sin(theta(1,i)); %Defining the sine of the joint angle

Rx_alpha=int64([1 0 0 0;0 ca -sa 0;0 sa ca 0;0 0 0 1]);  %Matrixof link twist
Dx_a=[1 0 0 a(1,i);0 1 0 0;0 0 1 0;0 0 0 1];             %Matrix of link length
Rz_theta=[ct -st 0 0;st ct 0 0;0 0 1 0;0 0 0 1];         %Matrixof joint angle
Dz_d=[1 0 0 0;0 1 0 0;0 0 1 d(1,i);0 0 0 1];             %Matrix of Joint offset

T{1,i}=simplify(Rz_theta*Dz_d*Dx_a*Rx_alpha);            %Constructing the Homogeneous transformation
end
end

