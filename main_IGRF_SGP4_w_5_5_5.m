
clc
clear all
close all

global Hs1 Hs2
Hs1 = 0;
Hs2 = 0;


global I ts mu0 start_date 
I = [2.22 0 0; 0 2.18 0; 0 0 0.5]*(10^(-2));

q = [1 0 0 0]';
w = [0.01;0.01;0.01]*pi/180;

ts = 0.1;
mu0 = 4*pi*10^(-7);

%-------------------------------------------------------------------------%
%                    TLE file name for running SGP4                       %
%-------------------------------------------------------------------------%
fname = 'tle_csswe.txt';
%simulation start date:  datenum(2012,09,14,00,59,48): 735126.041527778
start_date = datenum(2012,09,14,00,59,48);

time_total = 100 ;  

timespan = linspace(1,time_total,time_total);
t = 0;
for t_time = 0:0.1:time_total
     
    t = t + 1
    tsince = t_time/60;                  % tsince in minutes
    [pos,vel] = SGP4_eci( tsince, fname );   
    
%-------------------------------------------------------------------------%
%                        IGRF model of earth                           %
%-------------------------------------------------------------------------%    
  
    jd = datevec(addtodate(start_date,t_time,'second'));
    lat = eci2lla((pos*1000)',jd); % lat(1) = latitude, lat(2) = longitude, lat(3) = altitude

    B_eci = IGRF_eci(lat, jd);

    dcm = quat2dcm(q'); % convert quaternion to directional cosine matrix in body frame
    eul = quat2eul(q'); % convert quaternion to euler angles
   
    B_body = dcm*B_eci;

%-------------------------------------------------------------------------%
%                 residual magnetic moment torque                         %
%-------------------------------------------------------------------------%

%    m_res = [ 0.0059; 0.0083; -0.0004];
   T_res = hyst_trq2(B_body)
%     
% -------------------------------------------------------------------------%
%                       gravity gradient torque                            %
% -------------------------------------------------------------------------%
% 
Tgg = ComputeGGTorque(pos,q');

%-------------------------------------------------------------------------%
%                         aero dynamic torque                             %
%-------------------------------------------------------------------------%

Ta = ComputeAeroTorque(vel,q');
    

%-------------------------------------------------------------------------%
%                         calculating alpha angle                          %
%-------------------------------------------------------------------------%

%alpha = acosd(dot([1;0;0],S_body)/(norm(S_body)));
    
%-------------------------------------------------------------------------%
%                         Kinematic and Dynamic                           %
%-------------------------------------------------------------------------%
    
   T = Tgg + Ta + T_res;

    
    [q,w] = dynamic_kinematic_eq(w,q,T);
    
    w_x(t) = w(1);
    w_y(t) = w(2);
    w_z(t) = w(3);

    q1(t) = q(1);
    q2(t) = q(2);
    q3(t) = q(3);
    q4(t) = q(4);
    
    eul_1(t) = eul(1); 
    eul_2(t) = eul(2);
    eul_3(t) = eul(3);
    
    B_i_x(t) = B_eci(1);
    B_i_y(t) = B_eci(2);
    B_i_z(t) = B_eci(3);
    
    B_body_x(t) = B_body(1);
    B_body_y(t) = B_body(2);
    B_body_z(t) = B_body(3);
    
    Tx_res(t) = T_res(1);
    Ty_res(t) = T_res(2);
    Tz_res(t) = T_res(3);
    
    Tx_gg(t) = Tgg(1);
    Ty_gg(t) = Tgg(2);
    Tz_gg(t) = Tgg(3);
    
    Tx_a(t) = Ta(1);
    Ty_a(t) = Ta(2);
    Tz_a(t) = Ta(3);
    
    Tx(t) = T(1);
    Ty(t) = T(2);
    Tz(t) = T(3);

end
t = time_total*10 + 1;
s = 20;
timespan = linspace(1,t,t);
% time_x_label = store_t/6035;            %in days
time_x_label = timespan/(24*60*60*10);
figure(1)
plot(time_x_label, w_x*180/pi,'b',time_x_label,w_y*180/pi,'r', time_x_label,w_z*180/pi,'g' );
grid on;
title('angular rates');   xlabel('time (days)');  ylabel('angular rates (deg/sec)');
legend('x','y','z');
xlim([0 15]);
set(gca,'fontsize',s);

figure(2)
subplot(2,2,1);  plot(time_x_label,q1');   title('q0');  xlabel('time (days)'); ylabel('quaternion');
grid on;
set(gca,'fontsize',s);
subplot(2,2,2);  plot(time_x_label,q2');   title('q1');  xlabel('time (days)'); ylabel('quaternion');
grid on;
set(gca,'fontsize',s);
subplot(2,2,3);  plot(time_x_label,q3');   title('q2');  xlabel('time (days)'); ylabel('quaternion');
grid on;
set(gca,'fontsize',s);
subplot(2,2,4);  plot(time_x_label,q4');   title('q3');  xlabel('time (days)'); ylabel('quaternion');
grid on;
set(gca,'fontsize',s);

figure(3) 
plot(time_x_label, B_body_x,'b',time_x_label, B_body_y,'r',time_x_label, B_body_z,'g');
grid on;
title('Magnetic fields in body frame');   xlabel('time (days)');  ylabel('B (T)');
legend('x','y','z');
set(gca,'fontsize',s);

figure(4)
plot(time_x_label, B_i_x,'b',time_x_label, B_i_y,'r',time_x_label, B_i_z,'g');
grid on;
title('Magnetic fields in eci frame');   xlabel('time (days)');  ylabel('B (T)');
legend('x','y','z');
set(gca,'fontsize',s);


figure(5)
plot(time_x_label, Tx','b', time_x_label, Ty','r', time_x_label, Tz','g');  
grid on;
title('Total Torque');     xlabel('time (days)'); ylabel('Torque (Nm)');   
legend('x','y','z');
set(gca,'fontsize',s);

figure(6)
plot(time_x_label, Tx_gg','b', time_x_label, Ty_gg','r', time_x_label, Tz_gg','g'); 
title('Torque due to gravity gradient');    xlabel('time'); ylabel('Torque (Nm)');   
legend('x','y','z');
set(gca,'fontsize',s);

figure(7)
plot(time_x_label, Tx_a','b', time_x_label, Ty_a','r', time_x_label, Tz_a','g'); 
title('Torque due to aerodynamic drag');     xlabel('time'); ylabel('Torque (Nm)'); 
legend('x','y','z');
set(gca,'fontsize',s);

figure(8)
plot(time_x_label, Tx_res','b', time_x_label, Ty_res','r', time_x_label, Tz_res','g'); 
grid on;
title('Torque due to residual magnetic moment');     xlabel('time'); ylabel('Torque (Nm)'); 
legend('x','y','z');
set(gca,'fontsize',s);

figure(9)
plot(time_x_label, eul_1,'b',time_x_label, eul_2,'r',time_x_label, eul_3,'g');
grid on;
title('euler angles');       xlabel('time (days)'); ylabel('euler angles (rad)'); 
legend('x','y','z');
set(gca,'fontsize',s);

%figure(10)
%plot(time_x_label, alpha_plot);
%grid on;
%title('alpha plot');       xlabel('time (days)'); ylabel('alpha (deg)'); 
%set(gca,'fontsize',s);
