function [ wd ] = dyn_eq( w,T )

global I

%-------------------------------------------------------------------------%
%                    Dynamic equation                                     %
%-------------------------------------------------------------------------%

I_inv = inv(I);
I_w = I*w;                %Iw
cr_I_W = cross(w,I_w);   %Iw x w

wd = I_inv*(-cr_I_W + T);

end

