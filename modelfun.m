function dq = fcn1(t, q)
%nozzle velocity of bottle rocket

P_0 = 172370;       %Pa, initial internal pressure
P_2 = 101325;       %Pa, atmospheric pressure

D_b = 0.0592;       %m, diameter of bottle
D_n = 0.0218;       %m, diameter of nozzle
D_a = 0.0103;       %m, diameter of air sting
L_a = 0.197;        %m, length of air sting
L_s = 0.057;        %m, length of seal
A_b = pi*D_b^2/4;   %cm^2, area of bottle
A_n0 = pi*D_n^2/4;  %cm^2, area of nozzle
L_b = 0.2;          %m, length of bottle
L_0 = 0.065;         %m, initial height of water
L_2 = 0.0203;       %m, length of nozzle
V_0 = A_b*(L_b-L_0);%m^3, initial volume of air
rho_h2o = 1000;     %kg/m^3, density of water
rho_air = 1.225;    %kg/m^3, density of air (uncompressed)
mr = 0.06478;         %kg, mass of rocket without water
theta_0 = 60;         %degrees, launch angle from ground
g = 9.807;            %m/s^2, acceleration due to gravity
Cd = 0.3;           %, drag coefficient of rocket

v1 = q(1);
L1 = q(2);
vrx = q(3);
vrz = q(4);
x = q(5);
z = q(6);

lin = sqrt(x^2+z^2);    %linear distance traveled
vlin = sqrt(vrx^2+vrz^2);

if lin < L_a
    A_n = A_n0 - pi* D_a^2 / 4;
else
    A_n = A_n0;
end


v2 = A_b/A_n * v1;
P_1 =  (P_0 )* (V_0^1.4) / ((L_b-L1) * A_b)^1.4;

if(L1 <= 0) %bottle is ~empty
    dv1 = 0;
    dL1 = 0;
    L_2 = 0;
    v2 = 0;
else        %bottle contains water
    dL1 = - v1;
    dv1 = ((P_1)+rho_h2o/2*(v1^2-v2^2))/(rho_h2o*L1 + rho_h2o*A_b/A_n*L_2);
end

if (lin < L_s)%still on seal
    theta = theta_0;
    dvrx = cosd(theta)*(P_1) * A_n0 / (mr+rho_h2o*A_b*L1 + 0*rho_h2o*A_n*L_2);
    dvrz = sind(theta)*(P_1) * A_n0 / (mr+rho_h2o*A_b*L1 + 0*rho_h2o*A_n*L_2);
    dv1 = A_n / A_b * sqrt(dvrx^2 + dvrz^2);
else   %off of seal
    theta = atand(vrz/vrx);
    dv2 = A_b/A_n * dv1;
    dvrx = (cosd(theta)*-rho_h2o*A_b*(dL1*(vlin-v1)-L1*dv1) + cosd(theta)*rho_h2o*A_n*L_2*dv2 - cosd(theta)*rho_h2o*v2*A_n*(vlin-v2) - cosd(theta)*(Cd/2)*rho_air*A_b*vlin*abs(vlin))/(mr+rho_h2o*A_b*L1 + rho_h2o*A_n*L_2);
    dvrz = (sind(theta)*-rho_h2o*A_b*(dL1*(vlin-v1)-L1*dv1) + sind(theta)*rho_h2o*A_n*L_2*dv2 - sind(theta)*rho_h2o*v2*A_n*(vlin-v2) - (mr+rho_h2o*A_b*L1)*g - sind(theta)*(Cd/2)*rho_air*A_b*vlin*abs(vlin))/(mr+rho_h2o*A_b*L1 + rho_h2o*A_n*L_2);
 
end



if(z < -0.5)
    vrx = 0;
    vrz = 0;
    dvrz = 0;
    dvrx = 0;
end

dx = vrx;
dz = vrz;
dq = zeros(length(q),1);

dq(1) = dv1;
dq(2) = dL1;
dq(3) = dvrx;
dq(4) = dvrz;
dq(5) = dx;
dq(6) = dz;