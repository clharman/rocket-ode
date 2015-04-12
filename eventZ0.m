function [value,isterminal,direction] = eventZ0(t,y)
% Locate the time when the rocket passes the sting
value = y(6)+0.4;     % launch zero is greater than zero
isterminal = 1;   % Stop the integration
direction = 0;