function L = hyst_trq2(B_body)


len = 0.095;
dia = 0.001;
mu = (4*pi)*1e-7;

Vh = 6*pi*((dia/2)^2)*len;
Bs = 0.3;
Hc = 0.3381;
Br = 6.0618*1e-4;
% Bs = 0.74;
% Hc = 0.96;
% Br = 0.35;

global Hs1_x 
global Hs1_y
% Hc = 12;%constant characterising coercivity value of the hysteresis loop for each of the rod(A/m)
% Br = 0.004;
% %Hr = 0.004/(4*pi*1.5e-3);%constant characterising remanance value of the hysteresis loop for each of the rod(A/m)
% Bm = 0.025;%constant saturation value characterising the hysteresis major loop for each of the hysteresis rods(A/m)
% Vh = pi*15e-8;%volume of the hysteresis rod(m^3)
% 
% p = tan(pi*Br/(2*Bm))/Hc;
% 
% q0 = q(1,1);
% q1 = q(2,1);
% q2 = q(3,1);
% q3 = q(4,1);
% mu = (4*pi)*1e-7;
% Rbi = [q0^2+q1^2-q2^2-q3^2 2*(q1*q2+q0*q3) 2*(q1*q3-q0*q2); 2*(q1*q2-q0*q3) q0^2-q1^2+q2^2-q3^2 ...
%     2*(q2*q3+q0*q1); 2*(q1*q3+q0*q2) 2*(q2*q3-q0*q1) q0^2-q1^2-q2^2+q3^2];
% %Rib = Rbi';
% B = mag_field(t);
% Bb = Rbi*B;
% global Hs1 Hs2
% % H = d*Bb/mu;
%number_rods = [3 3 0];

mu_dash = (1.5*10^(4))/(1+(1.5*10^(4))/(((4*len)/(dia*sqrt(pi)))+2));
H_x = 3*B_body(1)/(mu);
H_y = 3*B_body(2)/(mu);

c1 = H_x-Hs1_x;
Hs1_x = H_x;

c2 = H_y-Hs1_y;
Hs1_y = H_y;


p = tan((pi*Br)/(2*Bs))/Hc;

if(c1<0)
    Bhyst_x = 2*Bs*atan(p*(H_x+Hc))/pi;
else 
    Bhyst_x = 2*Bs*atan(p*(H_x-Hc))/pi;
    
end

if(c2<0)
    Bhyst_y = 2*Bs*atan(p*(H_y+Hc))/pi;
else 
    Bhyst_y = 2*Bs*atan(p*(H_y-Hc))/pi;
    
end

Bhyst = [Bhyst_x, Bhyst_y, 0];
mhyst = Bhyst*Vh/(2*mu); % along each axis
L = cross(mhyst,B_body);
% 
% 
% d = [0 1 0];
% H = d*Bb/mu;
% c = H - Hs2;
% Hs2 = H;
% if(c<0)
%     Bhyst = 2*Bm*atan(p*(H+Hc))/pi;
% else 
%     Bhyst = 2*Bm*atan(p*(H-Hc))/pi;
%     
% end
% mhyst = Bhyst*Vh/mu;
% L = L1+ cross(mhyst*d',Bb);
end