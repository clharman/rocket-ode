[t,y] = ode45(@valfun,[0, 0.05], [0, 0.0793, 0, 0, 0, 0, 0]);

% time (s)
sampleT = [ 0.0000 0.0024 0.0048 0.0072 0.0096 0.0120 0.0144 ...
          0.0168 0.0192 0.0216 0.0240 0.0264 0.0288 0.0312 ...
          0.0336 0.0360 0.0384 0.0408 ];
% distance (m)
sampleX = [ 0.0000 0.0051 0.0102 0.0178 0.0279 0.0432 0.0584 ...
          0.0813 0.1092 0.1422 0.1803 0.2210 0.2667 0.3175 ...
          0.3658 0.4191 0.4750 0.5334 ];
% time step (s)
% fps = 1250; 3 is a magic number, taken from the spreadsheet
dT = 3 / 1250;
% calculate velocity
sampleV = zeros(1, length(sampleX)-1);
for i = 2 : length(sampleX) - 1
    % using central difference
    sampleV(i) = ( sampleX(i+1) - sampleX(i-1) ) / ( 2*dT ); 
end

linear_pos = sqrt(y(:,5).^2 + y(:,6).^2);
linear_vel = sqrt(y(:,3).^2 + y(:,4).^2);

%figure(1); plot(t,linear_pos); hold on; title('position'); plot(sampleT,sampleX, '*r');

%figure(2); plot(t,linear_vel); hold on; title('velocity'); plot(sampleT(1:end-1),sampleV, '*r');

%figure(3); plot(t,y(:,6)); title('height');

%figure(4); plot(t,y(:,5)); title('horizontal distance');

figure(5); plot(linear_pos,linear_vel); hold on; plot(sampleX(1:end-1),sampleV, '*r'); title('velocity vs position');