function [state_est, covar, measurements] = myParticleFilter(timeArray, x0, noiseModel, nParts, dt)

% Set seed
rng(100);

% Pre allocate data
weights = ones(nParts, 1)/nParts;

resampleThreshold = 1;

initialStatePDF = x0;

% Initial pdf - p(x) - I think this is it - past this should be p(x_k+1|x_k)

%Sample nParts particles from pdf of x - with replacement 
parts = datasample(pdf_x, nParts, 'Replacement', true);

for i = 1:timeArray
    % Do filtering
    
    
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


for i = 1:timeArray
    
    %Need to know p(y_k|x_k^i) fxn so that we can plug in particle i
    %generated
    pdf_yk_xki; 
    
    [posteriorDist, oldParticles, weights] = sampleFromProposalDist(parts,x0,noiseModel,pdf_yk_xki);
    
    % Renormalize weights
    weights = weights / sum(weights);
    
    %C: I commented this part out because according to table 3.4, this is not
    %done for bootstrap particle filter since it resamples every time 
    %{
    % Calculate effective # of particles (TODO: insert algo)
    Neff = 0;
    
    if Neff < resampleThreshold
        parts = bootStrapParticles(parts);
    end
    %}
    
    [parts, weights] = bootStrapParticles(parts, weights, 'myAlgo');
    
    % MMSE or MAP state estimate
    state_est = getEstimate(parts, linspace(min(parts), max(parts), nParts*100), 'MAP');
    
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
function [newParticles, newWeights] = bootStrapParticles(oldParticles, oldWeights, resamplingType)

switch resamplingType
    
case 'myAlgo'
        
    %Allocate space for resampled parameters
    newParticles = zeros(length(oldParticles), 1);

    % Sort weights of particles in descending order (first = max)
    [weight_cdf, oldWeightIndices] = sort(oldWeights, 'descend');

    for samples = 1:length(oldParticles)

        % Draw radom sample of 0 to 1 from uniform distribution
        u = rand;

        % Find first index where cdf is greater than Unif[0,1]
        part_index = find( weight_cdf > u , 1);

        %Add this to new particle list 
        newParticles(samples) = oldParticles(part_index);

    end

    % Uniform new weights
    newWeights = ones(length(oldWeight), 1) / length(oldWeights);

case 'stock'

    [newParticles, newWeights] = resample(oldParticles, oldWeights);

otherwise 
    
    error("Enter a valid resampling type to function: bootStrapParticles")
    
end

end
function [estimate, covar] = getEstimate(particles, support, typeOfEstimator)

switch typeOfEstimator
    
    case 'MMSE'
    
    %(TODO: insert algo)
    
    estimate = 1;
    covar    = 1;
    
    case 'MAP'
    
    % Kernel density estimate to generate a PDF
    density = ksdensity(particles, support, 'kernel', 'epanechnikov');
    
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
