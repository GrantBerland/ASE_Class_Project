function [observation_or_model] = measurement_model(state, mode)

if strcmp(mode, 'full')
    
    % Unpack position and convert to lat/lon/alt coordinates
    x = state(1);
    y = state(2);
    z = state(3);
    [latitude, longitude, altitude] = ecef2geod(x, y, z, 1e-3);
    
    % Get B-field state at current position
    time = datenum([2019 7 17 6 30 0]);
    B = igrf(time, latitude, longitude, altitude, 'geocentric');
    
    observation_or_model = state(7:9);
    
elseif strcmp(mode, 'linearized')

    observation_or_model = [zeros(3,8) ones(3,1)];

else
    error("Enter a valid mode to fnc:measurement_model!");
end
end