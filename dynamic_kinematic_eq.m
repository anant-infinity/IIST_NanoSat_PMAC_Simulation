function [ q,w ] = dynamic_kinematic_eq( w,q,T )

global I ts

%-------------------------------------------------------------------------%
%               Updating the angular rate w                               %
%-------------------------------------------------------------------------%

K1 = dyn_eq(w, T);
K2 = dyn_eq(w + 0.5*K1*ts, T);
K3 = dyn_eq(w + 0.5*K2*ts, T);
K4 = dyn_eq(w + K3*ts, T);

w_new = w + ((1/6)*ts*(K1 + 2*K2 + 2*K3 + K4));
w = w_new;

%-------------------------------------------------------------------------%
%               Updating the quaternion q                                 %
%-------------------------------------------------------------------------%

K1 = kine_eqn(q, w);
K2 = kine_eqn(q + 0.5*K1*ts, w);
K3 = kine_eqn(q + 0.5*K2*ts, w);
K4 = kine_eqn(q + K3*ts, w);

q_new = q + ((1/6)*ts*(K1 + 2*K2 + 2*K3 + K4));
q = q_new;

q = q/norm(q);                          %Normalization

if(q(1) < 0)
    q = -q;
end


end

