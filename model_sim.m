%This script runs modelfun with a set of initial conditions and
%  displays selected simulation results.
%Created for MECHENG 495 W14 Lab 4 by:
%  Section 003 Team 4: Colin Harman, Brian Freeburg, Joe Hendrickson

%CHANGE THESE-----------------------------------------------
L_0 = 0.064;            %m, initial height of water in bottle
theta_0 = 30;           %degrees, launch angle
psi = 30;

%ROCKET-SPECIFIC CONSTANTS----------------------------------
D_b = 0.0592;           %m, diameter of bottle
A_b = pi*D_b^2/4;       %cm^2, area of bottle
L_b = 0.196;            %m, length of bottle
V_0 = A_b*(L_b-L_0);    %m^3, initial volume of air in bottle

IC = [0, 0.064, 0, 0, 0, 0, psi*6894, theta_0*pi/180, V_0, 0, 0, 0];

%RUNS UNTIL ROCKET LEAVES STING
%pressure drop when leaving sting
options = odeset('Events',@eventLA);
[t1,y1] = ode45(@modelfun,[0, 6], IC,options);
y1(end,7) = y1(end,7) * 0.98;

%RUNS UNTIL ROCKET HITS GROUND
options = odeset('Events',@eventZ0,'InitialStep',t1(end)-t1(end-1));
[t2,y2] = ode45(@modelfun,[t1(end), t1(end)+30], y1(end,:), options);

%Concatenate arrays
t = cat(1, t1, t2);
y = cat(1, y1, y2);
    
linear_pos = sqrt(y(:,5).^2 + y(:,6).^2);
linear_vel = sqrt(y(:,3).^2 + y(:,4).^2);

%figure(1); plot(t,y(:,8));             title('angle vs t');
%figure(2); plot(t,linear_vel);         title('speed vs t');
%figure(3); plot(t,y(:,6));             title('height vs t');
%figure(3); plot(t,y(:,12));            title('y vs t');
%figure(4); plot(t,y(:,5));             title('horizontal distance vs t');
%figure(5); plot(y(:,5), y(:,6));       title('z vs x');
%figure(6); plot(t,y(:,3));             title('xvel vs time');
%figure(7); plot(t,y(:,4));             title('zvel vs time');
%figure(13); plot(t,y(:,11));           title('yvel vs time');
%figure(8); plot(t, y(:,8));            title('trajectory angle');
figure(10); plot(t, y(:,10)./(2*linear_vel*tand(9.6)/D_b));          title('spin rate');

%creates multiplot of x-z position, speed, angle, spin rate:
%{
figure % create new figure

subplot(1,4,1) % first subplot
plot(y(:,5), y(:,6))
%title('First subplot')
xlabel('x-position, m')
ylabel('z-position, m')
xlim([0 65])
ylim([-0.5 15])

subplot(1,4,2)
plot(t, y(:,8))
%title('First subplot')
ylabel('trajectory angle, rad')
xlabel('time, s')
%ylim([0 160])

subplot(1,4,3)
plot(t, linear_vel)
%title('First subplot')
ylabel('rocket speed, m/s')
xlabel('time, s')

subplot(1,4,4)
plot(t, y(:,10))
%title('First subplot')
ylabel('spin rate, rad/s')
xlabel('time, s')
ylim([0 160])
%}

y(end,5) * 3.28
