function [ B_eci ] = IGRF_eci( lat, jd )

% global start_date

% jd = datevec(addtodate(start_date,t,'minute'));
% 
% lat = eci2lla((pos*1000)',jd);

B_ned = igrf(datenum(jd),lat(1),lat(2),lat(3)/1000);

[B_ecef(1),B_ecef(2),B_ecef(3)] = ned2ecefv(B_ned(1),B_ned(2),B_ned(3),lat(1),lat(2));

%selecting the earth rotation angular velocity and deriving the dcm

DCM = dcmeci2ecef('IAU-2000/2006',jd);

% w = 7.2921159e-5;   %rad/sec
% time_num = addtodate(start_date,t,'minute');
% 
% R = [cos(w*(time_num)) sin(w*(time_num)) 0;
%     -sin(w*(start_date+t)) cos(w*(start_date+t)) 0;
%               0                       0          1];

B_eci=(DCM'*(B_ecef'))*10^(-9);

end

