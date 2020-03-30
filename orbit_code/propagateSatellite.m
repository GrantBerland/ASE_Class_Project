function [r,v] = propagateSatellite(r,v,dt)
% Unpack state variables
t0 = 0;
x0 = [r; v];
               
% Integrate EOMs, propagate satellite forward dt
[t, x] = ode45(@(t,x)orbitFnc(t,x), [t0 t0+dt], x0);

% Repack state variables to return
r = x(end,1:3);
v = x(end,4:6);
end
function dxdt = orbitFnc(t, x)
% x state contains:
% x(1:3) - position vector
% x(4:6) - velocity vector

% Earth gravitational parameter
mu = 3.986004418e5; % km^3/s^2

% Discrete time update under central body inverse r^2 force
dxdt    = zeros(6,1);
dxdt(1) = x(4); 
dxdt(2) = x(5); 
dxdt(3) = x(6); 
dxdt(4) = -mu/norm(x(1:3))^3*x(1);
dxdt(5) = -mu/norm(x(1:3))^3*x(2);
dxdt(6) = -mu/norm(x(1:3))^3*x(3);
end