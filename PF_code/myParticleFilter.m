function [state_est, covar, measurements] = myParticleFilter(timeArray, nParts, x0)

weights = zeros(nParts, 1);
parts   = zeros(nParts, 1);

resampleThreshold = 1;

initialState = x0;

% Generate initial particles (TODO: insert algo)
parts = ones(nParts, 1);

for i = 1:timeArray
    
    [posteriorDist, weights] = sampleFromProposalDist(parts, x0);
    
    % Renormalize weights
    weights = weights / sum(weights);
    
    % Calculate effective # of particles (TODO: insert algo)
    Neff = 0;
    
    if Neff < resampleThreshold
        parts = bootStrapParticles(parts);
    end
    
    % MMSE or MAP state estimate
    state_est = getEstimate(parts, 'MMSE');
    
end


end
function newParticles = bootStrapParticles(oldParticles, distribution)
%(TODO: insert algo)
end
function [posterior, w] = sampleFromProposalDist()
%(TODO: insert algo)
end
function [estimate, covar] = getEstimate(particles, typeOfEstimator)
%(TODO: insert algo)
if strcmp(typeOfEstimator, 'MMSE')
    estimate = 1;
    covar    = 1;
elseif strcmp(typeOfEstimator, 'MAP')
    estimate = 1;
    covar    = 1;
end

end
