[t,y] = ode45(@modelfun,[0, 30], [0, 0.065, 0, 0, 0, 0]);


linear_pos = sqrt(y(:,5).^2 + y(:,6).^2);
linear_vel = sqrt(y(:,3).^2 + y(:,4).^2);

%figure(1); plot(t,linear_pos); title('position');

%figure(2); plot(t,linear_vel); title('velocity vs t');

%figure(3); plot(t,y(:,6)); title('height vs t');

%figure(4); plot(t,y(:,5)); title('horizontal distance vs t');

figure(5); plot(y(:,5), y(:,6)); title('z vs x');

%figure(6); plot(t,y(:,3)); title('xvel vs time');
%figure(7); plot(t,y(:,4)); title('zvel vs time');

%figure(8); plot(t, atand(y(:,4)/y(:,3))); title('trajectory angle');

y(end,5)