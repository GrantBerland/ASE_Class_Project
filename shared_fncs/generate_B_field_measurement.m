function B = generate_B_field_measurement(noiseModel, pertMagnitude, plotOn, r, v, dt, measurement_length)

nSteps = measurement_length;

time = datenum([2019 7 17 6 30 0]);
tol = 1e-3;

timeStepFrac = 4;

% Position array over +/- 0.5 dt
x_array = linspace(r(1) - v(1)*dt/timeStepFrac, r(1) + v(1)*dt/timeStepFrac, nSteps);
y_array = linspace(r(2) - v(2)*dt/timeStepFrac, r(2) + v(2)*dt/timeStepFrac, nSteps);
z_array = linspace(r(3) - v(3)*dt/timeStepFrac, r(3) + v(3)*dt/timeStepFrac, nSteps);


[lat, lon, alt] = ecef2geod(x_array, y_array, z_array, tol);

% Gets nSteps points +/- 0.1 degrees ahead and behind spacecraft for 
% measurement model inversion


% Call to IGRF for a series of locations
B = igrf(time, lat, lon, alt, 'geocentric');

% Add perturbations according to noiseModel and pertMagnitude
%noiseModel    = 'gaussian';
%noiseModel    = 'gmm';
%noiseModel    = 'exp';
%noiseModel    = 'students-t';
%noiseModel    = 'none';

if strcmp(noiseModel, 'gaussian')
    pert = randn(3,nSteps)*pertMagnitude;
elseif strcmp(noiseModel, 'students-t')
    pert = trnd(3,[3 nSteps])*pertMagnitude;
elseif strcmp(noiseModel, 'gmm')
    mu = randn(5, 1)*10;
    gm = gmdistribution(mu,diag(pertMagnitude));
    pert = reshape(gm.random(3*nSteps)*pertMagnitude, 3, nSteps);
elseif strcmp(noiseModel, 'none')
    pert =  zeros(3,nSteps);
elseif strcmp(noiseModel, 'exp')
    % Two-sided exponential noise
    pert = -sign(randn(3,nSteps))*pertMagnitude.*log(rand(3,nSteps));
end

% Bz component typically ~10x larger in magnitude
B = B + pert'.*[1, 1, 10];

% Plot B-field components
if plotOn == 1
    figure(1); subplot(3,1,1); title('B-field Measurements');
    plot(B(:,1),'.'); ylabel('B_x [nT]'); title('True B-field'); grid on;
    subplot(3,1,2); plot(B(:,2),'.'); ylabel('B_y [nT]'); grid on;
    subplot(3,1,3); plot(B(:,3),'.'); ylabel('B_z [nT]'); grid on; xlabel('Sample');
end

end