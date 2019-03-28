% Function to Plot the ABB Robot Configuration based on the Joint angles
% The inputs are the numeric Joint Angles, T: Homogeneous transformations cell and
% HT is the total homogenous transformation.
% The output is the plot of the Robot configuration and composite HT
% This function is made for Course Project RBE501 
% made by Mostafa Atalla

function HT = plotarm(thetav,av,T,HT,x,y,index)
global theta a dx dy

DOF=size(thetav,2); %Compute the DOF based on the DH vectors size

for i=1:DOF+1
    T{index,i}=subs(T{index,i},[theta,a,dx,dy],[thetav,av,x,y]);
end
HT=double(subs(HT,[theta,a,dx,dy],[thetav,av,x,y]));

%% Plotting the Robot Configuration %%
% Initiation 
Or=zeros(3,DOF+2);                     % DH Frames origins Location Matrix
tipframe=zeros(4,4);                  % Tip Frame matrix includes the direction cosines of this frame

% Plot Settings
plot3(0,0,0);                          % Plot the origin in 3D Plot
daspect([1 1 1]);                      % Aspect Ratio 1
hold on
grid on
xlabel('X');
ylabel('Y');
zlabel('Z');

% Plot Base Frame 
plot3([0,20],[0,0],[0,0],'r','linewidth',1)     % Plot Base Frame X Axis
%text(20,0,0,'Xo')
plot3([0,0],[0,20],[0,0],'g','linewidth',1)     % Plot Base Frame Y Axis
%text(0,20,0,'Yo')
plot3([0,0],[0,0],[0,20],'b','linewidth',1)     % Plot Base Frame Z Axis
%text(0,0,20,'Zo')
plot3(0,0,0,'.','markersize',20)                % Plot Joint 1 

% Calculating the DH frames origins location
for i=1:DOF+1
    T1=T{index,i};
    Or(:,i+1)=T1(1:3,4);
    plot3([Or(1,i),Or(1,i+1)],[Or(2,i),Or(2,i+1)],[Or(3,i),Or(3,i+1)],'linewidth',4,'color','black');  %Plot each Link
    plot3(Or(1,i+1),Or(2,i+1),Or(3,i+1),'.','markersize',30)                                           %Plot the Joint Location  
end

% Calculating the tip Frame coordinate directions
for ii=1:3
    tipBasisVector=zeros(4,1);
    tipBasisVector(4,1)=1;
    tipBasisVector(ii,1)=15;
    tipframe(:,ii)=HT*tipBasisVector;
end

% Plotting the tip frame 
plot3([Or(1,DOF+2),tipframe(1,1)],[Or(2,DOF+2),tipframe(2,1)],[Or(3,DOF+2),tipframe(3,1)],'r','linewidth',1)
%text(tipframe(1,1),tipframe(2,1),tipframe(3,1),'Xtip')

plot3([Or(1,DOF+2),tipframe(1,2)],[Or(2,DOF+2),tipframe(2,2)],[Or(3,DOF+2),tipframe(3,2)],'g','linewidth',1)
%text(tipframe(1,2),tipframe(2,2),tipframe(3,2),'Ytip')

plot3([Or(1,DOF+2),tipframe(1,3)],[Or(2,DOF+2),tipframe(2,3)],[Or(3,DOF+2),tipframe(3,3)],'b','linewidth',1)
%text(tipframe(1,3),tipframe(2,3),tipframe(3,3),'Ztip')


end

