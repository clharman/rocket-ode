%have to change L_0 as well as IC L
%dist = [];
%j = 1;
L_b = 0.247;
psi = 150;
D_b = 0.1016;



%for i = 20:1:46
    L_0 = 0.112;
    theta_0 = 35;
    A_b = pi*D_b^2/4;
    V_0 = A_b*(L_b-L_0);
    IC = [0, L_0, 0, 0, 0, 0, psi*6894, theta_0*pi/180, V_0];

    
    options = odeset('Events',@eventZ0);
    [t,y] = ode45(@OCEfun,[0 30], IC, options);
%    dist(j, :) = [y(end,5)*3.28, i];
%    j = j+1;
%    y = [];
%end

%plot(dist(:,1),dist(:,2));title('angle (deg) vs total distance')

%figure(1); plot(t,linear_pos); title('position');

%figure(2); plot(t,linear_vel); title('velocity vs t');

%figure(3); plot(t,y(:,8)); title('angle vs t');

%figure(4); plot(t,y(:,5)); title('horizontal distance vs t');

figure(5); plot(y(:,5), y(:,6)); title('z vs x');

%figure(6); plot(t,y(:,3)); title('xvel vs time');
%figure(7); plot(t,y(:,4)); title('zvel vs time');

%figure(8); plot(t, atand(y(:,4)/y(:,3))); title('trajectory angle');
%figure(9); plot(t, y(:,7));

y(end,5) * 3.28
