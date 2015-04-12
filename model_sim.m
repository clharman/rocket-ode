%have to change L_0 as well as IC L

psi = 30;
IC = [0, 0.069, 0, 0, 0, 0, psi*6894];
m_r = 0.06478;
phi_fin = 9.6;
D_b = 0.0592;
J = 0.000000008;

%pressure drop when leaving sting
options = odeset('Events',@eventLA);
[t1,y1] = ode45(@modelfun,[0, 6], IC,options);
y1(end,7) = y1(end,7) * 0.95;
y1(end,1) = y1(end,1) * 0.5;

%convert linear to angular energy when water is gone
options = odeset('Events',@eventL1);
[t2,y2] = ode45(@modelfun,[t1(end), t1(end)+30], y1(end,:), options);
v2_0 = sqrt(y2(end,3)^2 + y2(end,4)^2);
theta = atan(y2(end,4)/y2(end,3));
v2_1 = sqrt(m_r * (v2_0^2) / (m_r + J/ (0.5*D_b*tand(phi_fin))^2 ));
%y2(end,3) = v2_1 * cos(theta);
%y2(end,4) = v2_1 * sin(theta);

options = odeset('Events',@eventZ0);
[t3,y3] = ode45(@modelfun,[t2(end), t2(end)+30], y2(end,:), options);


t = cat(1, t1, t2, t3);
y = cat(1, y1, y2, y3);


linear_pos = sqrt(y(:,5).^2 + y(:,6).^2);
linear_vel = sqrt(y(:,3).^2 + y(:,4).^2);

%figure(1); plot(t,linear_pos); title('position');

figure(2); plot(t,linear_vel); title('velocity vs t');

%figure(3); plot(t,y(:,6)); title('height vs t');

%figure(4); plot(t,y(:,5)); title('horizontal distance vs t');

figure(5); plot(y(:,5), y(:,6)); title('z vs x');

%figure(6); plot(t,y(:,3)); title('xvel vs time');
%figure(7); plot(t,y(:,4)); title('zvel vs time');

%figure(8); plot(t, atand(y(:,4)/y(:,3))); title('trajectory angle');
%figure(9); plot(t, y(:,7));

y(end,5) * 3.28
