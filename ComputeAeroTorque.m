function Ta = ComputeAeroTorque(v_eci, q)
%% Corrected Aerodynamic Torques - SMAD J Wertz

km2m    = 1e3;    
rho     = 1e-13;                                % Atmospheric density at an altitude of 700 km (in kg/m^3)
Cd      = 2.4;                                 % Aerodynamic drag coefficient (Typically between 2 and 2.5)
A       = 0.1*0.4;                              % Spacecraft area exposed to the atmosphere in the direction of the orbital velocity(in m^2)

V           = norm(v_eci)*km2m;                 % Norm of the spacecraft velocity in m/sec

dcm_body    =  quat2dcm(q);
v_body      =  dcm_body*v_eci;                  % Spacecraft velocity in the body frame of reference
uv          =  v_body/norm(v_body);             % Unit vector in the velocity direction

r =[2.601; -0.218; -8.086]*10^(-3);                    % Taken from CSSWE Thesis - Offset of the aerodynamic C.o.P from the C.o.M of the spacecraft

constant = 0.5*rho*(V^2)*Cd*A;

Ta = constant*cross(uv,r);                      % Computation of aerodynamic torque in N-m

end