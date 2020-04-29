function [state_est, covar, measurements] = myParticleFilter(timeArray,noiseModel, nParts, x0)

%set seed
rng(100);
%Pre allocate data
weights = zeros(nParts, 1);
parts   = zeros(nParts, 1);

resampleThreshold = 1;

initialState = x0;

%Initial pdf - p(x) - I think this is it - past this should be p(x_k+1|x_k)


%Sample nParts particles from pdf of x - with replacement 
parts = datasample(pdf_x,nParts);

% Generate initial particles (TODO: insert algo)
parts = ones(nParts, 1);


%Initial weights will be equal
w_part = ones(nParts,1)/nParts;

function [posterior, oldParticles, w] = sampleFromProposalDist(parts,x0,noiseModel,pdf_yk_xki)
%(TODO: insert algo)

for i = 1:nParts
    
%To sample for particles, first need to sample noise 
%This depends on what noise model we are using 

if strcmp(noiseModel, 'gaussian')
    noise_i = randn;
elseif strcmp(noiseModel, 'students-t')
    pert = trnd(3);
elseif strcmp(noiseModel, 'gmm')
    %C: This confuses me 
    mu = randn(5, 1)*10;
    gm = gmdistribution(mu,diag(pertMagnitude));
elseif strcmp(noiseModel, 'none')
    pert =  0;
elseif strcmp(noiseModel, 'exp')
    % Two-sided exponential noise
    % C: Also do not undestand
    pert = -sign(randn(3,nSteps))*pertMagnitude.*log(rand(3,nSteps));
end

%To sample can propagate particles k-1 and noise through dynamics
%C: I will leave this to you - need to save all particles - note we are
%doing one particle at a time 

%To find weights - use p(y_k|x_k^i) and evaluate at x_k^i

w(i) = pdf_yk_xki;

end
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
    
    parts = bootStrapParticles(oldParticles,weights);
    
    % MMSE or MAP state estimate
    state_est = getEstimate(parts, 'MMSE');
    
end


end

function [posterior, oldParticles, w] = sampleFromProposalDist(parts,x0,noiseModel,pdf_yk_xki)
%(TODO: insert algo)

for i = 1:nParts
    
%To sample for particles, first need to sample noise 
%This depends on what noise model we are using 

if strcmp(noiseModel, 'gaussian')
    noise_i = randn;
elseif strcmp(noiseModel, 'students-t')
    pert = trnd(3);
elseif strcmp(noiseModel, 'gmm')
    %C: This confuses me 
    mu = randn(5, 1)*10;
    gm = gmdistribution(mu,diag(pertMagnitude));
elseif strcmp(noiseModel, 'none')
    pert =  0;
elseif strcmp(noiseModel, 'exp')
    % Two-sided exponential noise
    % C: Also do not undestand
    pert = -sign(randn(3,nSteps))*pertMagnitude.*log(rand(3,nSteps));
end

%To sample can propagate particles k-1 and noise through dynamics
%C: I will leave this to you - need to save all particles - note we are
%doing one particle at a time 

%To find weights - use p(y_k|x_k^i) and evaluate at x_k^i

w(i) = pdf_yk_xki;

end
end
function newParticles = bootStrapParticles(oldParticles, oldWeights)
%(TODO: insert algo)

%Allocate space 
newParticles = ones(nParts,1);

%Make a CDF of oldWeights - need a value of CDF for each oldWeight ie
%particle - make sure the CDF is in order for each weight - need to make
%sure correct weight matches particle

cdf_weights;

for i = 1:nParts
    
    %Draw radom sample of 0 to 1 from uniform distribution
    rand_x = rand;
    
    %Find first index where cdf is greater than rand_x
    index_particle = find(cdf_weights>rand,1);
    
    %Add this to new particle list 
    newParticles(i) = oldParticles(index_particle);
end 

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
