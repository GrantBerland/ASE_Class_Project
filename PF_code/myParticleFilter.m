function [state_est, covar, measurements] = myParticleFilter(timeArray, x0, Q, noiseModel, pertMagnitude, nParts, dt)

state_est = zeros(length(x0), length(timeArray)); % 9 x tSpan
covar     = zeros(length(x0), length(x0), length(timeArray)); % 9 x 9 x tSpan

% Assign initial uniform weight to all particles
weights = ones(nParts, 1)/nParts;

resampleThreshold = 0.6;

% Initial pdf - p(x) - I think this is it - past this should be p(x_k+1|x_k)
initialPDF = normpdf(x0, Q);


%Sample nParts particles from pdf of x - with replacement 
parts = datasample(pdf_x, nParts, 'Replacement', true);

for i = 1:timeArray
    
    % Get posterior PDF 
    [posterior, oldParticles, w] = sampleFromProposalDist(parts, weights, x0, noiseModel, pdf_yk_xki);
    
    % Generate measurement
    B = generate_B_field_measurement(noiseModel, pertMagnitude, 0, r, v, dt, dt*100);
    
    % Get weights
    weights = computWeights();
    
    % Normalize weights s.t. sum(weights) == 1
    weights = weights / sum(weights);
    
    % Test if particles need to be resampled
    Neff = 1/sum(weights.^2);
    if Neff < resampleThreshold*nParts
        [parts, weights] = resample(parts, weights);
    end
    
    % Get estimate at current timestep
    [state_est(:,i), covar(:,:,i)] = getEstimate(parts, weights, support, 'MAP');
end

end
function [posterior, oldParticles, w] = sampleFromProposalDist(parts,x0,noiseModel,pdf_yk_xki)
%(TODO: insert algo)

for i = 1:nParts
    
%To sample for particles, first need to sample noise 
%This depends on what noise model we are using 

%To sample can propagate particles k-1 and noise through dynamics
%C: I will leave this to you - need to save all particles - note we are
%doing one particle at a time 

%To find weights - use p(y_k|x_k^i) and evaluate at x_k^i

w(i) = pdf_yk_xki;

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