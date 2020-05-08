

addpath('./EKF_code','./igrf_code','./orbit_code', './shared_fncs');
set(0,'DefaultFigureWindowStyle','docked')


% Establish satellite orbit (ISS orbit) and convert to (r0, v0)
[OE] = TLE2OrbitalElements(getSatelliteTLE(25544));
[satellite_r0,satellite_v0] = orb2rv(OE.a_km*(1-OE.e^2), OE.e, ...
               deg2rad(OE.i_deg),deg2rad(OE.Omega_deg),...
               deg2rad(OE.omega_deg),0,...
               0,0,0);
           
% Convert propagated position to lat/lon/alt coordinates
latLonAltTolerance = 1e-3;
[latitude, longitude, altitude] = ecef2geod(satellite_r0(1)*1000, ... % [deg, deg, m]
                                            satellite_r0(2)*1000,...
                                            satellite_r0(3)*1000, ...
                                            latLonAltTolerance);
% Get B-field at r0
time = datenum([2019 7 17 6 30 0]);
B0 = igrf(time, latitude, longitude, altitude, 'geodetic')';

% Noise models to add to IGRF
%noiseModel    = 'gaussian';   % AWGN
noiseModel    = 'gmm';        % Gaussian mixture model 
%noiseModel    = 'exp';        % exponential noise
%noiseModel    = 'students-t'; % student's-t noise
%noiseModel    = 'none';       % no noise

plotOn  = 1;
nStates = 9;

simLength = 1000;

dt = 10;
nParts = 100;
x0 = [satellite_r0, satellite_v0, B0];

[state_est, covar, meas] = myParticleFilter(linspace(0, simLength, simLength), x0, nParts, dt);


t_end = length(measurements);
t_span = linspace(0, t_end, length(state_est));


j = 3;
figure(j);
for i = 1:3   
  subplot(3,1,i); grid on; hold on;
  plot(t_span, state_est(:,i+3*(j-1)))
  ylabel(sprintf('B_{%i} [nT]', i));
end
subplot(3,1,1); title('Spacecraft B-field Measurements');

