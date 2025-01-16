function [value,isterminal,direction] = wavespeed_events_fourvar(~,u,~)
nx=100;
value       = u(1:end)-0.9;
isterminal  = [0*u(1:end-1);0.99];
direction   = 0*u;

end