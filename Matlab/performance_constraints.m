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
Md = diag([2*ones(1,3) 0.1*ones(1,3)]');  
Cd = diag([20*ones(1,3) 0.5*ones(1,3)]');  %Cartesian admittance parameters

%External force to be applied thoughout the simulation in XYZ directions (simulating human force)
F_h =[0 10 0 0 0 0]'; 
%F_h =[5 5 5]';  %alternative scenario
% F_h =[5 0 0]';  %alternative scenario
% F_h =[0 0 -5]';  %alternative scenario

%%%Performance constraints parameters%%%
separate_TR = 1; %0: use translational index for both, 1: use seperate indices
option = 2;      %0: no limits, 1: maniplty, 2: MSV, 3: inv. Condition Number
if (option<=1)      %maniplty
    w_crit = [0.01  0.2];
    w_th   = [0.025 0.5];
    lambda = [1   1];
elseif (option==2)  %MSV
    w_crit = [0.02  0.1];
    w_th   = [0.2   0.4];
    lambda = [1   1];
elseif (option==3)  %inv. Condition Number
    w_crit = [0.05  0.05];
    w_th   = [0.4  0.4];
    lambda = [1   1];
end

%initialisation of other parameters
sim_time = 4; %Simulation time in seconds
time = 0:Ts:sim_time;
q = q0';
A=zeros(6,1);
x_dotprev=zeros(6,1);
Q = []; W=[]; F_v=[]; P=[];
figure(1); subplot(3,2,[1 3 5]); 
r.plot(q0,'delay',Ts,'floorlevel',-0.4,'perspective')  %Create figure to simulate the robot motion in real time
% xlim([-1 1]); ylim([-1 1]); zlim([-0.4 1])

clc; disp('Calculating the simulation...')
for i=1:length(time)
    J       = r.jacob0(q); JT = J(1:3,:); JR = J(4:6,:);
    JT_aug  = JT*(eye(length(q0),length(q0))-pinv(JR)*JR); %Eq. 13
    JR_aug  = JR*(eye(length(q0),length(q0))-pinv(JT)*JT);
    
    Wall(i)     = sqrt(det(J*J'));
    Wt(i)       = sqrt(det(JT*JT')); %Translational manipulability in the weak sense
    Wr(i)       = sqrt(det(JR*JR'));
    Wt_aug(i)   = sqrt(det( JT_aug*JT' )); %TTM in the strong sense
    Wr_aug(i)   = sqrt(det( JR_aug*JR' ));
    [U,sigma,V] = svd(J);      S(:,i) = (diag(sigma));
        detUV(:,i)= [det(U); det(V) ];
    [U,sig_t,V] = svd(JT_aug); St(i)  = min(svd(JT)); %MSV translational J only
        detUVt(:,i) = [det(U); det(V) ]; 
    [U,sig_r,V] = svd(JR_aug); Sr(i)  = min(svd(JR)); %MSV rotational J only
        detUVr(:,i) = [det(U); det(V) ]; 
    St_aug(i)   = min(diag(sig_t)); %MSV augmented translational J
    Sr_aug(i)   = min(diag(sig_r)); %MSV augmented rot J
    CN(i)       = min(diag(sigma))/(max(diag(sigma)));
    CNt(i)      = min(diag(sig_t))/max(diag(sig_t));
    CNr(i)      = min(diag(sig_r))/max(diag(sig_r));
    

    %Select criterion
    switch(option)
        case {0, 1} %manipulability
            if (separate_TR);	W(:,i)= [Wt_aug(i) Wr_aug(i)]';
            else                W(i)  = Wt_aug(i);        end
        case 2  %MSV
            if (separate_TR);	W(:,i)= [St_aug(i) Sr_aug(i)]';
            else                W(i)  = min(S(:,i));        end
        case 3  %Condition number
            if (separate_TR);	W(:,i)= [CNt(i) CNr(i)]';
            else                W(i)  = CN(i);        end
    end

    %for translation (or both when separate_TR=0)
    if ( W(1,i) > w_th(1) )
       k(1,i)=0;
    else
       k(1,i)=lambda(1)*(1/(W(1,i)-w_crit(1))-1/(w_th(1)-w_crit(1)));
    end
        
    if (separate_TR) %calculate additional rotations
        if ( W(2,i) > w_th(2) )
           k(2,i)=0;
        else
           k(2,i)=lambda(2)*(1/(W(2,i)-w_crit(2))-1/(w_th(2)-w_crit(2)));
        end
    end
    
    A = -grad_w(q,r,J,option,separate_TR);

    if (separate_TR)
        F_v(1:3,i) = k(1,i)*A(1:3); %Eq.(1)
        F_v(4:6,i) = k(2,i)*A(4:6); %Eq.(1)
    else
        F_v(:,i) = k(i)*A; %Eq.(1)
    end

    if (option==0)
        F_v(:,i) = zeros(6,1);
    end

    %Cartesian admittance
    x_dot(:,i)= (Md/Ts + Cd) \ (Md*x_dotprev/Ts + F_h - F_v(:,i) ) ;
    Jinv = pinv(J);
  
    q_dot = Jinv*(x_dot(:,i)); %differential inverse kinematics 

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
    x_dotprev=x_dot(:,i); %save for next iteration
    Q_dot(:,i)=q_dot;  %save joint reference velocities
    
    %show the progress
    if (mod(i,0.05/Ts) ==0) 
       clc; disp([num2str(round(100*i/length(time))) '%' ])
    end
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
scatter3(P(1,:),P(2,:),P(3,:),10,W(1,:))
colormap(Colormp/max(max(Colormp)))
% colorbar('location','EastOutside'); %Sometimes this messes with the plot. Try at your own risk
caxis([0,max(W(1,:))])