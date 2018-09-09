%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Gravity-Gradient Torques as a function of the
% current orbital position.
% Input             - r_eci    - Satellite's current orbital position vector in km
%                   - q        - Satellite's current attitude quaternion
% Output            - Tgg      - Gravity Gradient Torque in N-m
% Global Variable   - I        - Spacecraft Inertia matrix
% Local Variables   - km2m     - Conversion of kilometres to metres
%                   - Mu       - Earth's Gravitational constant parameter
%                   - R        - Magnitude(norm) of the position vector(r_eci) in metres
%                   - dcm_body - DCM corresponding to the current spacecraft attitude quaternion ( ECI to the body frame )
%                   - r_body   - Components of the satellite's current orbital position vector expressed in body frame in km
%                   - R0       - Normalised, unit r_body vector
% Function Call     - my_const, q2dcm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tgg = ComputeGGTorque(r_eci,q)

global  I

%%Etika's Code
km2m        =  1000;
Mu          =  6.67259*10^(-11) * 5.98*10^(24);


R           =  norm(r_eci)*km2m;
dcm_body    =  quat2dcm(q);
r_body      =  dcm_body*r_eci;
R0          =  r_body/norm(r_body);
Tgg         =  (3*Mu/R^3)*cross(R0,I*R0);


end