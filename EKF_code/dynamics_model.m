function [state_or_model] = dynamics_model(state, mode, dt)

% Earth gravitational parameter
mu = 3.986004418e5; % km^3/s^2

if strcmp(mode, 'full')
    % Unpack spacecraft position-velocity state
    r = state(1:3);
    v = state(4:6);
    
    % Propagate satellite forward one time step
    [r, v] = propagateSatellite(r, v, dt);
    
    % Convert propagated position to lat/lon/alt coordinates
    latLonAltTolerance = 1e-3;
    [latitude, longitude, altitude] = ecef2geod(r(1), r(2), r(3), latLonAltTolerance);
    
    % Get B-field at next time step
    time = datenum([2019 7 17 6 30 0]);
    B = igrf(time, latitude, longitude, altitude, 'geodetic');
    
    % Return propagated state
    state_or_model = [r, v, B]';
    
elseif strcmp(mode, 'linearized')
    % Unpack spacecraft position-velocity-B-field state
    r = state(1:3);
    v = state(4:6);
    B = state(7:9);
  
    % Further unpack state position
    x = r(1); y = r(2); z = r(3);
    
    % Calculate dv/dr
    dvdr  = [-(mu*(- 2*x^2 + y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2), ...
              (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2),...
              (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);...
              (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2),...
              -(mu*(x^2 - 2*y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2),...
              (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);...
              (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2),...
              (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2),...
              -(mu*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^(5/2)];
    
    % Calculate dB/dr = -d^2 V/dr^2
    dBdr = [eye(3)];
          
    state_or_model = [eye(3), zeros(3), zeros(3);...
                      dvdr  , eye(3)   , zeros(3);...
                      dBdr  , zeros(3) , eye(3)];
                  
else
    error("Enter a valid mode to fnc:dynamics_model!");
end
end