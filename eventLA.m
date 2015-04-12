function [value,isterminal,direction] = eventLA(t,y)
% Locate the time when the rocket passes the sting
value = sqrt(y(5)^2+y(6)^2)-0.197;     % Detect sting passed
isterminal = 1;   % Stop the integration
direction = 0;