function dq = fcn1(t, q)

D_b = 0.0592;       %m, diameter of bottle
D_n = 0.0218;       %m, diameter of nozzle
D_a = 0.0103;       %m, diameter of air sting
L_a = 0.197;        %m, length of air sting
L_s = 0.057;        %m, length of seal
A_b = pi*D_b^2/4;   %cm^2, area of bottle
A_n0 = pi*D_n^2/4;  %cm^2, area of nozzle
L_b = 0.196;        %m, length of bottle
L_0 = 0.069;         %m, initial height of water
L_2 = 0.0203;       %m, length of nozzle

rho_h2o = 1000;     %kg/m^3, density of water
rho_air = 1.225;    %kg/m^3, density of air (uncompressed)
m_r = 0.06478;         %kg, mass of rocket without water
theta_0 = 30;         %degrees, launch angle from ground
g = 9.807;            %m/s^2, acceleration due to gravity
Cd_lock = 0.35;       %, drag coefficient of rocket when still
Cd_spin = 0.2;       %, drag coefficient of rocket when spinning

mu = 0.00089;       %Pa s, dynamic viscosity of water at 25 C
%eps_pl = 0.0000025; %mm, roughness of plastic
%rpr_1 = eps_pl/D_b; %-, relative pipe roughness, bottle
%rpr_2 = eps_pl/D_n; %-, relative pipe roughness, nozzle
K_minor = 0*0.08;    %-, resistance coefficient of bottle-nozzle contraction
phi_fin = 9.6;      %degrees, angle of fins

v_1 = q(1); %velocity of flow in the bottle
L_1 = q(2); %height of water in the bottle
vrx = q(3); %velocity of rocket in x
vrz = q(4); %velocity of rocket in z
x = q(5);   %rocket position in x
z = q(6);   %rocket position in z
P_0 = q(7); %initial gage pressure in pa

V_0 = A_b*(L_b-L_0);%m^3, initial volume of air


    
%change from still to spinning drag coefficient
if(L_1 > 0)    
    Cd = Cd_lock;
else
    Cd = Cd_spin;
end

lin = sqrt(x^2+z^2);    %linear distance traveled
vlin = sqrt(vrx^2+vrz^2);

if lin < L_a
    A_n = A_n0 - pi* D_a^2 / 4;
else
    A_n = A_n0;
end

v_2 = A_b/A_n * v_1;
P_1 =  (q(7))* (V_0^1.4) / ((L_b-L_1) * A_b)^1.4;

Re_1 = rho_h2o*v_1*D_b/mu;
Re_2 = rho_h2o*v_2*D_n/mu;
f_d_1 = 0*64/Re_1;
f_d_2 = 0*64/Re_2;

if(L_1 <= 0) %bottle is ~empty
    dv_1 = 0;
    dL_1 = 0;
    L_2 = 0;
    v_2 = 0;
else        %bottle contains water
    %change in height of water
    dL_1 = - v_1;
    %change in flow 
    %                                           major head loss in bottle
    dv_1 = ((P_1)+rho_h2o/2*(v_1^2-v_2^2) - rho_h2o*f_d_1*L_1/D_b*(v_1^2)/2 - rho_h2o*f_d_2*L_2/D_n*(v_2^2)/2 - rho_h2o*K_minor*(v_2^2)/2 ) / (rho_h2o*L_1 + rho_h2o*A_b/A_n*L_2);
end

if (lin < L_s)%still on seal
    theta = theta_0;
    dvrx = cosd(theta)*(P_1) * A_n0 / (m_r+rho_h2o*A_b*L_1 + rho_h2o*A_n*L_2);
    dvrz = sind(theta)*(P_1) * A_n0 / (m_r+rho_h2o*A_b*L_1 + rho_h2o*A_n*L_2);
    dv_1 = A_n / A_b * sqrt(dvrx^2 + dvrz^2);
else   %off of seal
    theta = atand(vrz/vrx);
    dv2 = A_b/A_n * dv_1;
    dvrx = (cosd(theta)*-rho_h2o*A_b*(dL_1*(vlin-v_1)-L_1*dv_1) + cosd(theta)*rho_h2o*A_n*L_2*dv2 - cosd(theta)*rho_h2o*v_2*A_n*(vlin-v_2) - cosd(theta)*(Cd/2)*rho_air*A_b*vlin*abs(vlin))/(m_r+rho_h2o*A_b*L_1 + rho_h2o*A_n*L_2);
    dvrz = (sind(theta)*-rho_h2o*A_b*(dL_1*(vlin-v_1)-L_1*dv_1) + sind(theta)*rho_h2o*A_n*L_2*dv2 - sind(theta)*rho_h2o*v_2*A_n*(vlin-v_2) - (m_r+rho_h2o*A_b*L_1)*g - sind(theta)*(Cd/2)*rho_air*A_b*vlin*abs(vlin))/(m_r+rho_h2o*A_b*L_1 + rho_h2o*A_n*L_2);
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

dq(1) = dv_1;
dq(2) = dL_1;
dq(3) = dvrx;
dq(4) = dvrz;
dq(5) = dx;
dq(6) = dz;
dq(7) = 0;