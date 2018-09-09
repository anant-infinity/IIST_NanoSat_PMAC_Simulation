function [ qd ] = kine_eqn( q, w )

%-------------------------------------------------------------------------%
%                    Kinematic equation                                   %
%-------------------------------------------------------------------------%

qk = [q(1) -q(2) -q(3) -q(4);
      q(2)  q(1) -q(4)  q(3);
      q(3)  q(4)  q(1) -q(2);
      q(4) -q(3)  q(2)  q(1)];

wk = [0;w(1);w(2);w(3)];

qd = 0.5*qk*wk;

end
