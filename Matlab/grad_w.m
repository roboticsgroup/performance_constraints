function A = grad_w(q,r,J,option,separate_TR)
%Approximate the gradient of the manipulability index with respect to the
%Cartesian frame of the end effector (ONLY in 1 direction. Faster approach)
%
%   Arguments:
%   q0:     current position
%   r:      robot handle
%   J:      Jacobian in q0 used for the control
%   option: Which criterion to use
%   separate: set to 1 when using separate translational and rotational indices
%
% UPDATE 06/2020. This is an improved and fasterunpublished version

dx=1e-5;    %infinitesimal step in positive direction
A = zeros(6,1);


for dir=1:6 %for each direction
    JT = J(1:3,:); JR = J(4:6,:); 
    JT_aug = JT*(eye(r.n,r.n)-pinv(JR)*JR);
    JR_aug = JR*(eye(r.n,r.n)-pinv(JT)*JT);

    %Make a temporary Jacobian...
    if (dir<=3) Jc = JT_aug; %...for translation
    else        Jc = JR_aug; %...for rotation
    end

    %calculate current performance metric
    switch(option)
        case {0, 1}      %manipulability
            if (separate_TR);   w_q0 = sqrt(det( Jc*Jc' ));
            else            w_q0 = sqrt(det( J*J' ));        end
        case 2  %MSV
            if (separate_TR);   w_q0 = min(svd(Jc));
            else            w_q0 = min(svd(J));        end
        case 3  %Condition number
            if (separate_TR);   sigma=svd(Jc);w_q0 = min(sigma)/(max(sigma));
            else            sigma=svd(J); w_q0 = min(sigma)/(max(sigma));        end
    end

    %calculate neighbour performance metric
    dp=[0 0 0 0 0 0]';
    dp(dir)=dx;
    qv=q + pinv(J)*dp; %Eq.3
    Jv = r.jacob0(qv);
    JT = Jv(1:3,:); JR = Jv(4:6,:);
    JT_aug = JT*(eye(r.n,r.n)-pinv(JR)*JR);
    JR_aug = JR*(eye(r.n,r.n)-pinv(JT)*JT);

    %Make a temporary Jacobian...
    if (dir<=3) Jc = JT_aug; %...for translation
    else        Jc = JR_aug; %...for rotation
    end

    switch(option)
        case {0, 1} %manipulability
            if (separate_TR);   w_next = sqrt(det( Jc*Jc' ));
            else            w_next = sqrt(det( Jv*Jv' )); end
        case 2  %MSV
            if (separate_TR);   w_next = min(svd(Jc));
            else            w_next = min(svd(Jv)); end
        case 3  %Condition number
            if (separate_TR);   sigma=svd(Jc);w_next = min(sigma)/(max(sigma));
            else            sigma=svd(Jv); w_next = min(sigma)/(max(sigma));        end
    end
    A(dir)=(w_next-w_q0) / dx;  %Eq.4
end

