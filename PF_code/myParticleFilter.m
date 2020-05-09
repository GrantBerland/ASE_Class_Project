function [state_est, covar, measurement_array] = myParticleFilter(timeArray, x0, Q, noiseModel, pertMagnitude, nParts, dt)

state_est = zeros(3, length(timeArray)); % 9 x tSpan
covar     = zeros(3, length(timeArray)); % 9 x 9 x tSpan
measurement_array = [];
% Assign initial uniform weight to all particles
weights = ones(nParts, 1)/nParts;

resampleThreshold = 0.7;

% Initialize particles normally about the B-field state
initParts = randn(nParts, 3)*10 + x0(7:end)';
incr = 0;
parts = initParts;
x = x0;
for i = 1:length(timeArray)
    
    r = x(1:3);
    v = x(4:6);
    
    %figure(1); hold on; grid on;
    %p1=scatter3(parts(:,1), parts(:,2), parts(:,3), weights*5000, '.');
    
    % Predict state with particle draw from importance PDF
    [parts, r, v] = sampleFromProposalDist(parts, weights, [r; v], noiseModel, pertMagnitude, dt);
    
    incr = incr + 1;
    if mod(incr, 10) == 0
        cla;
    end
    
    %p2=scatter3(parts(:,1), parts(:,2), parts(:,3), weights*5000, '.');
    xlabel('B_x'); ylabel('B_y'); zlabel('B_z');
    % Generate measurement
    measurement = Mag_Vec_Calibration_Carolina_shortcut_IGRF(noiseModel,pertMagnitude,0,[r; v],dt,dt*100);
    
    measurement(:,1) = measurement(:,1)+9e-6;
    measurement(:,2) = measurement(:,2)-11e-6;
    measurement(:,3) = measurement(:,3)+11.75e-5;
    %p3=scatter3(mean(measurement(:,1)), mean(measurement(:,2)), mean(measurement(:,3)), 'o');
    % Get weights based on observational ltikelihood PDF
    weights = computeWeights(mean(measurement, 1), parts, weights, r, v, dt);
    
    
    % Get and store estimate at current timestep
    [state_est(:,i), covar(:,i)] = getEstimate(parts, weights, 'MMSE');
    
    % Test if particles need to be resampled
    Neff = 1/sum(weights.^2);
    if Neff <= resampleThreshold * nParts
        % Draw new particles from previous weighted population
        parts = resample(parts, weights);
        
        % Reset weights uniformity
        %weights = ones(length(weights), 1)/length(weights);
    end
    
    %p4=scatter3(state_est(1,i), state_est(2,i), state_est(3,i), '*');
    %legend([p1(1),p2(1),p3(1),p4(1)],'x_{k-1}','x_{k}','Measurement','State Estimate')
    
    measurement_array = [measurement_array; mean(measurement, 1)];
end

end
function [newParts, r, v] = sampleFromProposalDist(parts, weights, x0_r_v, noiseModel, pertMagnitude, dt)

newParts = zeros(size(parts));
nSteps = length(parts);

% Noise models to sample from
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


for i = 1:length(parts)
    
    propagatedState = dynamics_model([x0_r_v; parts(i,:)'], 'full', dt);
    
    r = propagatedState(1:3);
    v = propagatedState(4:6);
    newParts(i,:) = propagatedState(7:9);
    
%To sample for particles, first need to sample noise 
%This depends on what noise model we are using 

%To sample can propagate particles k-1 and noise through dynamics
%C: I will leave this to you - need to save all particles - note we are
%doing one particle at a time 

%To find weights - use p(y_k|x_k^i) and evaluate at x_k^i
end

%Sample nParts particles from pdf of x - with replacement 
newParts = datasample(newParts + pert', length(parts), 1, 'Weights', weights, 'Replace', true);

end
function newWeights = computeWeights(measurement, parts, weights, r, v, dt)

innovation = measurement - parts;

newWeights = weights.*normpdf(mean(innovation, 2)*500000);

% Normalize weights s.t. sum(weights) == 1
newWeights = newWeights./sum(newWeights);

end
function [estimate, covar] = getEstimate(particles, weights, typeOfEstimator)

switch typeOfEstimator
    case 'MMSE'
        % Weighted average of all particles
        estimate = sum(weights.*particles);
        covar    = var(weights.*particles);
    case 'MAP'  
        estimate = [0;0;0];
        for i = 1:3
        % Kernel density estimate to generate a PDF
        density(i) = ksdensity(weights.*particles(:,i), ...
                            'kernel', 'epanechnikov');

        % Locate mode of posterior
        [~, idx] = max(density);

        % Estimate state and covariance within the 
        % resolution of the passed-in support
        estimate(i) = support(idx);
        covar    = var(density);
        end
    otherwise
        error("Enter a valid estimator type to function: getEstimate()")
end

end