function A = grad_w(q0,r,dir,J)
%Approximate the gradient of the manipulability index with respect to the
%Cartesian frame of the end effector (each direction seperately)
%
%   Arguments:
%   q0:     current position
%   r:      robot handle
%   dir:    direction to approximate
%   J:      Jacobian in q0 used for the control

v = 1;   %unit virtual velocity
Dt=0.01;    %integration interval
direc = [1 -1]; %positive and negative direction
n = length(q0);

w_q0 = sqrt(det( J*J' )); %Eq.15

for m = 1:2 %for each direction (positive and negative)
    x_dotv=[0 0 0]';
    x_dotv(dir)=v*direc(m);
    qv=q0 + pinv(J)*x_dotv*Dt; %Eq.3
    Jv = r.jacobn(qv);
    JTv = Jv(1:3,:); 
    JRV = Jv(4:6,:);
    w_next=sqrt(det( JTv*(eye(n,n)-pinv(JRV)*JRV)*JTv' ));
    
    if (w_next > w_q0)
        Delta_w(m)=w_next-w_q0;  %Eq.4
    else    %Exception handling Eq.6
        Delta_w(m)=0;
    end  

end

[A, index_of_sign] = max(Delta_w);
A = A*direc(index_of_sign); %Eq.5
end


 
 
