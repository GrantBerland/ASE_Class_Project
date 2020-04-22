function [observation_or_model] = measurement_model(state, mode, measurement_length, dt)

if strcmp(mode, 'full')
    
    timeStepFrac = 4;
    
    % Unpack position and convert to lat/lon/alt coordinates
    vx_arr = linspace(-state(4)/timeStepFrac, state(4)/timeStepFrac, measurement_length);
    vy_arr = linspace(-state(5)/timeStepFrac, state(5)/timeStepFrac, measurement_length);
    vz_arr = linspace(-state(6)/timeStepFrac, state(6)/timeStepFrac, measurement_length);
    
    x = state(1) + vx_arr * dt;
    y = state(2) + vy_arr * dt;
    z = state(3) + vz_arr * dt;
    
    % Vector (x,y,z)->(lat,long,alt) method
    [latitude, longitude, altitude] = ecef2geod(x, y, z, 1e-3);
    
    % Get B-field state at current position
    time = datenum([2019 7 17 6 30 0]);

    % Get B-field measurements over +/- one timestep
    observation_or_model = igrf(time, latitude, longitude, altitude, 'geocentric');
    
elseif strcmp(mode, 'linearized')

    % Directly observable Bx, By, Bz measurements from the state
    observation_or_model = repmat([zeros(3,6) eye(3)], measurement_length, 1);

else
    error("Enter a valid mode to fnc:measurement_model!");
end
end