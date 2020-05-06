function B = generate_B_field_dynamics(noiseModel, pertMagnitude, plotOn, lat, lon, alt)
nSteps = 1000;

% TODO: replace these with derived orbital parameters
% [x, y, z] = geod2ecef(latitude, longitude, altitude);
% [latitude, longitude, altitude] = ecef2geod(x, y, z, tol);
time = datenum([2019 7 17 6 30 0]);
lat  = linspace(-80, 80, nSteps);% deg
lon  = linspace(0, 120, nSteps); % deg
alt  = 500*1e3;                  % m

% Call to IGRF for a series of locations
[Bx, By, Bz] = igrf(time, lat, lon, alt, 'geocentric');

% Add perturbations according to noiseModel and pertMagnitude
%noiseModel    = 'gaussian';
%noiseModel    = 'gmm';
%noiseModel    = 'exp';
%noiseModel    = 'students-t';
%noiseModel    = 'none';

if strcmp(noiseModel, 'gaussian')
    pert = randn(3,nSteps)*pertMagnitude;
elseif strcmp(noiseModel, 'students-t')
    pert = trnd(3,nSteps)*pertMagnitude;
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

Bx = Bx + pert(1,:)';
By = By + pert(2,:)';
Bz = Bz + pert(3,:)';

% Plot B-field components
if plotOn == 1
    figure(1); subplot(3,1,1); plot(Bx); ylabel('B_x [nT]'); title('True B-field'); grid on;
    subplot(3,1,2); plot(By); ylabel('B_y [nT]'); grid on;
    subplot(3,1,3); plot(Bz); ylabel('B_z [nT]'); grid on; xlabel('Sample');
end
B = [Bx, By, Bz];
end