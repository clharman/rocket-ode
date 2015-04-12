function [value,isterminal,direction] = eventL1(t,y)
% Locate the time when height of water in rocket passes through zero 
value = y(2)-0;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;