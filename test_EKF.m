

addpath('./EKF_code','./igrf_code','./orbit_code');

% Establish satellite orbit (ISS orbit) and convert to (r0, v0)
[OE] = TLE2OrbitalElements(getSatelliteTLE(25544));
[satellite_r0,satellite_v0] = orb2rv(OE.a_km*(1-OE.e^2), OE.e, ...
               deg2rad(OE.i_deg),deg2rad(OE.Omega_deg),...
               deg2rad(OE.omega_deg),0,...
               0,0,0);
           
% Convert propagated position to lat/lon/alt coordinates
latLonAltTolerance = 1e-3;
% [deg, deg, m]
[latitude, longitude, altitude] = ecef2geod(satellite_r0(1)*1000, ... 
                                            satellite_r0(2)*1000, ...
                                            satellite_r0(3)*1000, ...
                                            latLonAltTolerance);

% Get B-field at r0
time = datenum([2019 7 17 6 30 0]);
B0 = igrf(time, latitude, longitude, altitude, 'geodetic')';

% Noise models to add to IGRF
%noiseModel    = 'gaussian';   % AWGN
%noiseModel    = 'gmm';        % Gaussian mixture model 
noiseModel    = 'exp';        % exponential noise
%noiseModel    = 'students-t'; % student's-t noise
%noiseModel    = 'none';       % no noise

plotOn  = 1;
nStates = 9;
nObservables = 3;
pertMagnitude = 5e-5;
%pertMagnitude = 0;

simulationLength = 200;
t_span = linspace(0, simulationLength, simulationLength);

Q = eye(nStates);

Q(7,7) = Q(9,9)*1e-6;
Q(8,8) = Q(9,9)*1e-6;
Q(9,9) = Q(9,9)*1e-7;

R = eye(nObservables)*2e-3;

R(3,3) = 0.0225;

dt = 1;

x0 = [satellite_r0; satellite_v0; B0];
[state_est, covar, measurements] = myEKF(t_span, Q, R, x0, dt, noiseModel, pertMagnitude);

residuals = measurements - state_est(:, 7:9);

figure(1); clf;
for i = 1:3   
  subplot(3,1,i); grid on; hold on;
  plot(t_span, measurements(:,i), 'b');
  plot(t_span, state_est(:,i+3*2), 'r');
  plot(t_span, state_est(:,i+3*2)+covar(:,i+3*2,i+3*2), '--g','Linewidth',0.5);
  plot(t_span, state_est(:,i+3*2)-covar(:,i+3*2,i+3*2), '--g','Linewidth',0.5)
  ylabel(sprintf('B_{%i} [nT]', i));
  if i == 1
      legend('True B-field','Estimated B-field','\pm 1 \sigma');
  end
end
subplot(3,1,1); title('Spacecraft B-field Measurements');

figure(2); clf;
for i = 1:3   
  subplot(3,1,i); grid on; hold on;
  %plot(t_span, residuals(:, i))
  histogram(residuals(:, i), 100);
  ylabel(sprintf('B_{%i} [nT]', i));
end
subplot(3,1,1); title('Spacecraft B-field Measurement Residuals');

%{
figure(); hold on;
[X,Y,Z] = sphere;
surf(X*6378,Y*6378,Z*6378);
scatter3(state_est(:,1), state_est(:,2), state_est(:,3), '.');
axis equal;
%}
