% Hand Exoskeleton - Course Project - RBE 501 %
% Made by Mostafa Atalla %

close all
clc
clear

%% Initialization %%
% Symbolic Parameters 
global theta a dx dy
syms t thetaCMC(t) thetaMCPAD(t) thetaMCPFE(t) thetaPIP(t) thetaDIP(t)     %Symbolic thetas
syms theta1dd(t) theta2dd(t) theta3dd(t) thetaMCPdot(t) thetaPIPdot(t) thetaDIPdot(t)
syms m1 m2 m3 I1 I2 I3
syms Lmc Lpp Lmp Ldp                                                       %Symbolic Bones Lengthes
syms dx dy                                                                 %Symbolic translational for carpal bones


% Symbolic DH Parameters Vectors and Rototranslational matrix
a= [Lmc 0 Lpp Lmp Ldp];                                                    %Link length
alpha= [-pi/2 pi/2 0 0 0];                                                 %Joint twist
d= [0 0 0 0 0];                                                            %Link twist
theta= [thetaCMC(t) thetaMCPAD(t) thetaMCPFE(t) thetaPIP(t) thetaDIP(t)];  %Joint Angle
Troto=[0 0 1 dx;1 0 0 dy;0 1 0 0;0 0 0 1];                                 %Rototranslational Matrix 
                                                                           %X coordinates of the center of mass of each link: Proximal, Middle and Distal

% Numeric Values for plotting
NOF=4;
DOF=size(theta,2);
thetaCMCn=[0 0 0 0];                                                       %Carpometacarpal joint angles vector from index to little
thetaMCPADn=[0 0 0 0];                                                     %Metacarpophalangeal Abduction/Adduction joint angles vector from index to little
thetaMCPFEn=[-40 -40 -40 -40];                                             %Metacarpophalangeal Flexion/Extension joint angles vector from index to little
thetaPIPn=[-40 -40 -40 -40];                                               %Proximal Interphalangeal joint angles vector from index to little
thetaDIPn=[-20 -20 -20 -20];                                               %Distal Interphalangeal joint angles vector from index to little

Lmcn=[68 64.6 58 53.69];                                                   %Metacarpal bones length vector from index to little
Lppn=[39.78 44.63 41.37 32.74];                                            %Proximal Interphalangeal Bones Length vector from index to little
Lmpn=[22.38 26.33 25.65 18.11];                                            %Middle Interphalangeal Bones Length vector from index to little
Ldpn=[15.82 17.4 17.3 15.96];                                              %Distal Interphalangeal Bones Length vector from index to little

x=[-25 0 25 50];                                                           %Fingers translational shift components in X Vector
y=[10 20 20 15];                                                           %Fingers translational shift components in Y Vector

% Numeric DH Paramters Vectors
theta_numeric=deg2rad([thetaCMCn;thetaMCPADn;thetaMCPFEn;thetaPIPn;thetaDIPn]);     %Numeric Joint Angles Matrix
a_numeric=[Lmcn;0 0 0 0;Lppn;Lmpn;Ldpn];                                            %Numeric Link Length Matrix

%% Kinematics: Compute the Transformation matrices w.r.t base frame and the total transformation
% T is the cell containing the transformation of each frame with respect to
% the base frame.
% HT is the total transformation from the base frame to the finger tip
T=cell(NOF,DOF+1);
HT=cell(1,NOF);
for i=1:NOF
    T{i,1}=Troto;
    CT= dhparam2matrix(theta,d,a,alpha);
    for ii=1:DOF
        T{i,ii+1}=CT{ii};
    end
    for ii=1:DOF
        T{i,ii+1}=simplify(T{i,ii}*T{i,ii+1});                             %Multiplying the H.Transformations to compute the totaltransformation
    end
    HT{i}=T{i,ii+1};
end

%% Kinematics: Compute the jacobian matrices applicable to each finger
%to be adjusted to cover the last three links (finger links)
for i=1:NOF
    J=jacobiancom(T,i);                                                    %The jacobian cell includes the jacobian matrix for each finger from index to little finger.
end

%% Plotting the hand configuration
figure(1);
for i=1:NOF
    thetav=theta_numeric(:,i)';
    av=a_numeric(:,i)';
    P = plotarm(thetav,av,T,HT{i},x(i),y(i),i);
title('Hand Configuration')
end

%% Dynamics: compute the transformations for the last three links
%Consider the 3 links of the fingers only and place the origin frame at the
%MCP joint and neglect the adduction/abduction motion.
% T is the cell containing the transformation of each frame with respect to
% the base frame.
% HT is the total transformation from the base frame to the finger tip
Tdyn=cell(NOF,3);
HTdyn=cell(1,NOF);
for i=1:NOF
    CT= dhparam2matrix(theta(1,3:5),d(1,3:5),a(1,3:5),alpha(1,3:5));
    for ii=1:3
        Tdyn{i,ii}=CT{ii};
    end
    for ii=1:2
        Tdyn{i,ii+1}=simplify(Tdyn{i,ii}*Tdyn{i,ii+1});                    %Multiplying the H.Transformations to compute the totaltransformation
    end
    HTdyn{i}=Tdyn{i,ii+1};
end

%% Dynamics: Jacobian 
%Obtaining the translational part of the kineatic energy
for i=1:NOF
    Jdyn=jacobiancom(Tdyn,i);                                              %The jacobian cell includes the jacobian matrix for each finger from index to little finger.
end

%% Calculating the Kinetic Energy
%Calculating the jacobian that corresponds to each center of mass
Jcom=cell(1,3);
CoM=[Lpp/2 0 0;Lpp/2 Lmp/2 0;Lpp/2 Lmp/2 Ldp/2];

q=[thetaMCPFE thetaPIP thetaDIP];
qdot=[thetaMCPdot thetaPIPdot thetaDIPdot];

for i=1:3
    Jcom{i}=subs(Jdyn(1:3,1:3),[Lpp Lmp Ldp],CoM(i,:));
end

% Calculating the total Kinetic Energy
KE_trans=0.5*qdot*(m1*Jcom{1}.'*Jcom{1}+m2*Jcom{2}.'*Jcom{2}+m3*Jcom{3}.'*Jcom{3})*qdot.';
KE_Rot=0.5*I1*thetaMCPdot^2+0.5*I2*(thetaMCPdot+thetaPIPdot)^2+0.5*I3*(thetaMCPdot+thetaPIPdot+thetaDIPdot)^2;
KE=KE_trans+KE_Rot;

% Computing the dynamic Equations using Lagrange
u=sym(zeros(3,1));
u(1,1)=diff((functionalDerivative(KE,thetaMCPdot)),t)-(functionalDerivative(KE,thetaMCPFE));
u(2,1)=diff((functionalDerivative(KE,thetaPIPdot)),t)-(functionalDerivative(KE,thetaPIP));
u(3,1)=diff((functionalDerivative(KE,thetaDIPdot)),t)-(functionalDerivative(KE,thetaDIP));
u=subs(u,[diff(thetaMCPdot(t), t) diff(thetaPIPdot(t), t) diff(thetaDIPdot(t), t) diff(thetaMCPFE(t), t) diff(thetaPIP(t), t) diff(thetaDIP(t), t)],[theta1dd theta2dd theta3dd thetaMCPdot thetaPIPdot thetaDIPdot]);
u=simplify(u);  

% Use collect function to separate the terms of the accelerations to be
% able to obtain the dynamics 