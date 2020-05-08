function [pert_state_est, covar, measurement_array] = myEKF(timeArray, Q, R, x0, dt, noiseModel, pertMagnitude)

% Number of states [x,y,z,vx,vy,vz,Bx,By,Bz]
n = 9;
nMeasurements = dt*100; % 100 Hz sampling rate

pert_state_est = zeros(length(timeArray), n);
covar          = zeros(length(timeArray), n, n);
measurement_array = [];

initial_state = x0;
initial_covar = eye(n)*1e0;

for i = 1:length(timeArray)
    
    % EKF prediction step
    [x_minus, P_minus] = EKF_predict(Q, initial_state, initial_covar, dt);
    
    % Simulated measurement at the spacecraft location
    measurement = generate_B_field_measurement(noiseModel, pertMagnitude, 0, ...
        initial_state(1:3), initial_state(4:6), dt, nMeasurements);
    initial_state(4:6)
    % EKF measurement update step
    [pert_state_est(i,:), covar(i,:,:)] = EKF_update(R, measurement, x_minus, P_minus, dt);
    
    % Set current step state estimate to be initial guess for next step
    initial_state = pert_state_est(i,:)';
    initial_covar = reshape(covar(i,:,:), 9, 9);
    
    measurement_array = [mean(measurement); measurement_array];
end
end
function [state_est, covar] = EKF_predict(Q, initial_state, initial_covar, dt)

% Predict state through nonlinear dynamics
state_est = dynamics_model(initial_state, 'full', dt);

% Linearize dynamics about previous measurement-updated state
F = dynamics_model(initial_state, 'linearized', dt);

% Estimate error covariance
covar     = F*initial_covar*F' + Q;

end
function [state_est, covar] = EKF_update(R, measurement, initial_state, initial_covar, dt)

% Compute linearized measurement model about current state
H = measurement_model(initial_state, 'linearized', length(measurement), dt);

% Repeat the R matrix according to the measurement size
R = diag(repmat(diag(R), length(measurement), 1));

% Compute Kalman gain
K_gain = initial_covar * H' / ( H * initial_covar * H' + R);

% Compute updated state measurement
innovation = measurement - measurement_model(initial_state, 'full', length(measurement), dt);
state_est = initial_state + K_gain * reshape(innovation', [], 1);

% Compute updated covariance
covar = (eye(length(state_est)) - K_gain * H) * initial_covar;

end