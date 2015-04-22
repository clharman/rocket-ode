 function dq = simplefun(t, q) %This function contains differential 
 %  equations that describe the behavior of a water bottle rocket.
 %This is a simplified version of modelfun.m
 %Created for MECHENG 495 W14 Lab 4 by:
 %  Section 003 Team 4: Colin Harman, Brian Freeburg, Joe Hendrickson
 %Written by Colin Harman with guidance from GSI Colin Jiang, 4/21/2015
 
 
%--------------------------------------------------------------------------
%-------------CONSTANTS----------------------------------------------------
%Rocket dimensions
L_b = 0.196;                %m, length of bottle
L_2 = 0.0203;               %m, length of nozzle
D_b = 0.0592;               %m, diameter of bottle
D_n = 0.0218;               %m, diameter of nozzle
A_b = pi*D_b^2/4;           %m^2, area of bottle (calculated from D_b)
A_n0 = pi*D_n^2/4;          %m^2, area of nozzle (calculated from D_n)

%Other rocket characteristics
m_r = 0.06478;              %kg, mass of rocket without water
Cd_spin = 0.21;             %-, drag coefficient of rocket when spinning
Cd_y = 1.17;                %-, drag coefficient for side of rocket

%Launcher dimensions
D_a = 0.0103;               %m, diameter of air sting                       
L_a = 0.209;                %m, length of air sting
L_s = 0.057;                %m, length of seal
A_s = pi*D_a^2/4;           %cm^2, area of air sting (calculated from D_a)

%Physical constants
rho_h2o = 997;             %kg/m^3, density of water
rho_air = 1.225;            %kg/m^3, density of air (uncompressed)
g = 9.807;                  %m/s^2, acceleration due to gravity

%Environmental variables
vy_wind = 0;             %m/s, wind speed in the y direction (2.23 m/s = 5 mph)
vx_wind = 0;                %m/s, wind speed in the x direction

%--------------------------------------------------------------------------
%-------------STATE VARIABLES----------------------------------------------
v_1 = q(1);                 %velocity of flow in the bottle
L_1 = q(2);                 %height of water in the bottle
vrx = q(3);                 %velocity of rocket in x
vrz = q(4);                 %velocity of rocket in z
x = q(5);                   %rocket position in x
z = q(6);                   %rocket position in z
P_0 = q(7);                 %initial gage pressure in pa (so that a pressure drop can be affected externally)
theta = q(8);               %pointing angle of the rocket in radians
V_0 = q(9);                 %initial volume of air (so that it does not need to be changed inside the function)
omega = q(10);              %spin speed of the rocket in rad/s
vry = q(11);                %velocity of rocket in y
y = q(12);                  %rocket position in y


%--------------------------------------------------------------------------
%-------------CONSTANT ALGEBRAIC EXPRESSIONS-(never change)----------------
lin = sqrt(x^2+z^2);        %magnitude of distance from launch point
vlin = sqrt(vrx^2+vrz^2);   %magnitude of x-z velocity
V = (L_b-L_1)*A_b;          %air volume
P = P_0*V_0^1.4/V^1.4;      %air pressure
asx = vrx - vx_wind;        %airspeed in x direction
asy = vry - vy_wind;        %airspeed in y direction
aslin = sqrt(asx^2+vrz^2);  %magnitude of x-z windspeed
Cd = Cd_spin;

%--------------------------------------------------------------------------
%-------------LAUNCH CONDITIONS--------------------------------------------
if lin < L_s                %rocket is on the launcher seal
    on_seal = true;         %   switch
else
    on_seal = false;
end

if lin < L_a                %rocket is on the air sting
    A_n = A_n0-A_s;         %   effective nozzle area is smaller
else
    A_n = A_n0;
end

if(L_1 <= 0)                %all fuel is gone
    empty = true;           %   switch
else
    empty = false;
end


%--------------------------------------------------------------------------
%-------------VARIABLE ALGEBRAIC EXPRESSIONS-(depend on launch conditions)-
%Constant flow rate
v_2 = A_b/A_n * v_1;


%--------------------------------------------------------------------------
%-------------DERIVATIVE DEFINITIONS---------------------------------------
if empty
    dv_1 = 0;
    dL_1 = 0;
    L_2 = 0;
    v_2 = 0;
else
    dL_1 = -v_1; 
    dv_1 = (P+rho_h2o/2*v_1^2*(1-(A_b/A_n)^2)) / (rho_h2o*L_1 + rho_h2o*A_b/A_n*L_2);
end

if on_seal
    dtheta = 0;
    dvrx = cos(theta)*P*A_n0 / (m_r+rho_h2o*A_b*L_1 + rho_h2o*A_n*L_2);
    dvrz = sin(theta)*P*A_n0 / (m_r+rho_h2o*A_b*L_1 + rho_h2o*A_n*L_2);
    dv_1 = A_n / A_b * sqrt(dvrx^2 + dvrz^2);
    dvry = 0;
else
    dv2 = A_b/A_n * dv_1;
    dvrx = (cos(theta)*-rho_h2o*A_b*(dL_1*(vlin-v_1)-L_1*dv_1) + cos(theta)*rho_h2o*A_n*L_2*dv2 - cos(theta)*rho_h2o*A_n*(vlin-v_2)*v_2 - cos(theta)*(Cd/2)*rho_air*A_b*aslin^2)/(m_r+rho_h2o*A_b*L_1 + rho_h2o*A_n*L_2);
    dvrz = (sin(theta)*-rho_h2o*A_b*(dL_1*(vlin-v_1)-L_1*dv_1) + sin(theta)*rho_h2o*A_n*L_2*dv2 - sin(theta)*rho_h2o*A_n*(vlin-v_2)*v_2 - (m_r+rho_h2o*A_b*L_1)*g - sin(theta)*(Cd/2)*rho_air*A_b*aslin^2)/(m_r+rho_h2o*A_b*L_1 + rho_h2o*A_n*L_2);
    dvry = 0.5*rho_air*Cd_y*(asy)^2*D_b*L_b/(m_r+rho_h2o*A_b*L_1 + rho_h2o*A_n*L_2);
    dtheta = (asx*dvrz - vrz*dvrx)/(asx^2 + vrz^2);
end

dx = vrx;
dz = vrz;
dy = vry;

%--------------------------------------------------------------------------
%-------------STATE DERIVATIVES--------------------------------------------
dq = zeros(length(q),1);

dq(1) = dv_1;
dq(2) = dL_1;
dq(3) = dvrx;
dq(4) = dvrz;
dq(5) = dx;
dq(6) = dz;
dq(7) = 0;              %initial pressure does not change internally
dq(8) = dtheta;
dq(9) = 0;              %initial volume does not change internally
dq(10) = domega;
dq(11) = dvry;
dq(12) = dy;