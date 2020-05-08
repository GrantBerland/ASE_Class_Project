function [state_est, covar, measurements] = myParticleFilter(timeArray, x0, Q, noiseModel, pertMagnitude, nParts, dt)

state_est = zeros(length(x0), length(timeArray)); % 9 x tSpan
covar     = zeros(length(x0), length(x0), length(timeArray)); % 9 x 9 x tSpan

% Assign initial uniform weight to all particles
weights = ones(nParts, 1)/nParts;

resampleThreshold = 0.6;

% Initialize particles normally about the B-field state
initParts = (randn(nParts, 3) + x0(7:end)')*10;

parts = initParts;
x = x0;
for i = 1:timeArray
    
    r = x(1:3);
    v = x(4:6);
    
    % Predict state with particle draw from importance PDF
    parts = sampleFromProposalDist(parts, weights, [r; v], noiseModel, pertMagnitude, dt);
    
    % Generate measurement
    measurement = generate_B_field_measurement(noiseModel, pertMagnitude, 0, r, v, dt, dt*100);
    
    % Get weights based on observational likelihood PDF
    weights = computeWeights(measurement, parts, weights, r, v, dt);
    
    % Normalize weights s.t. sum(weights) == 1
    weights = weights / sum(weights);
    
    % Get and store estimate at current timestep
    [state_est(:,i), covar(:,:,i)] = getEstimate(parts, weights, support, 'MAP');
    
    % Test if particles need to be resampled
    Neff = 1/sum(weights.^2);
    if Neff <= resampleThreshold*nParts
        % Draw new particles from previous weighted population
        parts = resample(parts, weights);
        
        % Reset weights uniformity
        weights = ones(length(weights), 1)/length(weights);
    end
    
end

end
function [newParts, x_r_v] = sampleFromProposalDist(parts, weights, x0_r_v, noiseModel, pertMagnitude, dt)
%(TODO: insert algo)

newParts = zeros(size(parts));
nSteps = 1;

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
    
    x_r_v      = propagatedState(1:6);
    newParts(i,:) = propagatedState(7:9);
    
%To sample for particles, first need to sample noise 
%This depends on what noise model we are using 

%To sample can propagate particles k-1 and noise through dynamics
%C: I will leave this to you - need to save all particles - note we are
%doing one particle at a time 

%To find weights - use p(y_k|x_k^i) and evaluate at x_k^i
end

%Sample nParts particles from pdf of x - with replacement 
newParts = datasample(newParts, length(parts), 1, 'Weights', weights, 'Replace', true);

end
function newWeights = computeWeights(measurement, parts, weights, r, v, dt)

newWeights = zeros(length(weights), 1);
for i = 1:length(weights)

    newWeights(i) = measurement_model([r; v], 'full', length(measurement), dt);
    
end


end
function [estimate, covar] = getEstimate(particles, weights, support, typeOfEstimator)

switch typeOfEstimator
    case 'MMSE'
        % Weighted average of all particles
        estimate = sum(weights.*particles);
        covar    = var(weights.*particles);
    case 'MAP'  
        % Kernel density estimate to generate a PDF
        density = ksdensity(weights.*particles, support, ...
                            'kernel', 'epanechnikov');

        % Locate mode of posterior
        [~, idx] = max(density);

        % Estimate state and covariance within the 
        % resolution of the passed-in support
        estimate = support(idx);
        covar    = var(density);
    otherwise
        error("Enter a valid estimator type to function: getEstimate()")
end

end