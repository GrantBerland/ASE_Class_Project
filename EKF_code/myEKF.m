function [pert_state_est, covar] = myEKF(measurements, Q, R, r0, v0, B0, dt)

% Number of states [x,y,z,vx,vy,vz,Bx,By,Bz]
n = 9;

pert_state_est = zeros(length(measurements), n);
covar     = zeros(length(measurements), n, n);

initial_state = [r0 ; v0; B0];
initial_covar = eye(n)*1e5;

for i = 1:length(measurements)
    
    % EKF prediction step
    [x_minus, P_minus] = EKF_predict(Q, initial_state, initial_covar, dt);
    
    % EKF measurement update step
    [pert_state_est(i,:), covar(i,:,:)] = EKF_update(R, measurements(i,:), x_minus, P_minus);
    
    % Set current step state estimate to be initial guess for next step
    initial_state = pert_state_est(i,:)';
    initial_covar = reshape(covar(i,:,:), 9, 9);
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
function [state_est, covar] = EKF_update(R, measurement, initial_state, initial_covar)

% Compute linearized measurement model about current state
H = measurement_model(initial_state, 'linearized');

% Compute Kalman gain
K_gain = initial_covar * H' / ( H * initial_covar * H' + R);

% Compute updated state measurement
state_est = initial_state + K_gain * (measurement' - measurement_model(initial_state, 'full'));

% Compute updated covariance
covar = (eye(length(state_est)) - K_gain * H) * initial_covar;

end