% "Manipulator Performance Constraints on Cartesian Admittance Control for
% human-robot Cooperation"
%
% Description:
% Using manipulability index to constrain the robot above a certain
% threshold.
%
% Usage:
% 1. Set the desired parameters for the robot and the performance
% constraints
% 2. Run the simulation
% 3. View the simulated motion of the robot and the plots
%
% Requirements:
% - Robotics Toolbox for Matlab (http://www.petercorke.com/Robotics_Toolbox.html)
% 
% Notes: 
% - The simulation does not consider the joint limits so the robot
% might behave weird
% - This code has been tested with Robotics Toolbox 9.10 and Matlab R2014a
%
% Authors: Fotios Dimeas, Charalambos Papakonstantinou
%
% Copyright 2015 Fotios Dimeas
%
%%
clear all;clc; clf;
load('colormp');

%%%  Robot Parameters  %%%
r = lwr();	%create robot object using Robotics Toolbox
q0=[0 0.3491 0 -1.3963 0 1.3963 0];	%set to an initial configuration
Ts=0.01; %Set simulation step (in seconds)
joint_velocity_limit = 2; %rad/s (from the specs of LWR)
Md = diag([2*ones(1,3)]');  Cd = diag([50*ones(1,3)]');  %Cartesian admittance parameters

%External force to be applied thoughout the simulation in XYZ directions (simulating human force)
F_h =[0 10 0]'; 
%F_h =[5 5 5]';  %alternative scenario
% F_h =[5 0 0]';  %alternative scenario
% F_h =[0 0 -5]';  %alternative scenario

%%%Performance constraints parameters%%%
enable_constraints = 1; %0: constraints are disabled, 1: constraints are enabled
sim_time = 4; %Simulation time in seconds
w_crit=0.01;
w_th=0.025;
lambda = 100;


%initialisation of other parameters
time = 0:Ts:sim_time;
q = q0';
A=zeros(3,1);
x_dotprev=zeros(3,1);
Q = []; W=[]; F_v=[]; P=[];
figure(1); subplot(3,2,[1 3 5]); 
r.plot(q0,'delay',Ts,'floorlevel',-0.4,'perspective')  %Create figure to simulate the robot motion in real time
% xlim([-1 1]); ylim([-1 1]); zlim([-0.4 1])

clc; disp('Calculating the simulation...')
for i=1:length(time)
    J = r.jacobn(q);
    JT = J(1:3,:); 
    JR = J(4:6,:);
    JT_aug = JT*(eye(length(q0),length(q0))-pinv(JR)*JR); %Eq. 13
    W = [W sqrt(det( JT_aug*JT_aug' ))]; %Translational manipulability in the strong sense

    % Eq.(2) asymptotic increase
    if ( W(:,i) > w_th )
       k=0;
    else
       k=lambda*(1/(W(i)-w_crit)-1/(w_th-w_crit));
       for dir=1:3 %for each of the translational directions
             A(dir)= grad_w(q,r,dir,JT_aug);  
       end 
    end

    if (enable_constraints)
        F_v(:,i) = k*A; %Eq.(1)
    else
        F_v(:,i) = zeros(3,1);
    end

    %Cartesian admittance
    x_dot=inv( Md/Ts + Cd)*(Md*x_dotprev/Ts + F_h + F_v(:,i) ) ;

    %differential inverse kinematics
    q_dot = pinv(JT_aug)*x_dot;
    
    %check the velocity limits
    if (max(abs(q_dot)) > joint_velocity_limit)
        disp('Joint velocity limits exceeded! Simulation stopped!')
        break
    end

    %integration
    q=Ts*q_dot+q;

    hom = r.fkine(q);
    P(:,i) = [hom(1,4) hom(2,4) hom(3,4)]';
    Q(:,i)=q; %save joint positions
    x_dotprev=x_dot; %save for next iteration
    Q_dot(:,i)=q_dot;  %save joint reference velocities
end
disp('Done!')

%% Plots
disp('Watch simulation in figure 1...')
figure(1)
subplot(3,2,[1 3 5])
% r.plot(q0,'delay',Ts,'floorlevel',-0.4,'perspective')
r.animate(Q')
 
subplot(3,2,2)
plot(time(1:length(Q_dot)),Q_dot)
xlabel('time')
ylabel('rad/s')
title('Joint Velocities')
legend('1','2','3','4','5','6','7')
ylim([-joint_velocity_limit joint_velocity_limit]); %LWR joint velocity limits

subplot(3,2,4)
plot(time(1:length(F_v)),F_v)
xlabel('time')
ylabel('N')
title('Constraining Force (XYZ)')
legend('x','y','z')

subplot(3,2,6)
plot(time(1:length(P)),P)
xlabel('time')
ylabel('m')
title('Cartesian position')
legend('x','y','z')

subplot(3,2,[1 3 5])
hold on
% r.plot(q0) %reset to initial position
scatter3(P(1,:),P(2,:),P(3,:),10,W)
colormap(Colormp/max(max(Colormp)))
% colorbar('location','EastOutside'); %Sometimes this messes with the plot. Try at your own risk
caxis([0,max(W)])