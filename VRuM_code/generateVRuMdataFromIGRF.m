function [B_tot_P1, B_ext_vec] = generateVRuMdataFromIGRF(nSteps, noiseModel, pertMagnitude)

fs       = 100; % Hz (100 Hz in Swarm doc)
duration = 10; % s

% Get IGRF solutions at nSteps points
time       = datenum([2019 7 17 6 30 0]); % datetime
latitude   = linspace(79, 80, nSteps);    % deg
longitude  = linspace(119, 120, nSteps);  % deg
altitude   = 500*1e3;                     % m

% External magnetic field measurement
B_external = igrf(time, latitude, longitude, altitude, 'geodetic')';

% Noise model to add to external B-field
if strcmp(noiseModel, 'gaussian')
    pert = randn(3,nSteps)*pertMagnitude;
elseif strcmp(noiseModel, 'students-t')
    pert = trnd(3,nSteps)*pertMagnitude;
elseif strcmp(noiseModel, 'gmm')
    mu = randn(5, 1)*10;
    gm = gmdistribution(mu,diag(pertMagnitude));
    pert = reshape(gm.random(3*nSteps)*pertMagnitude, 3, nSteps);
elseif strcmp(noiseModel, 'exp')
    % Two-sided exponential noise
    pert = -sign(randn(3,nSteps))*pertMagnitude.*log(rand(3,nSteps));
elseif strcmp(noiseModel, 'none')
    pert =  zeros(3,nSteps);
end

B_external = B_external + pert;

% Generate total B-field components
e_B = [ .8  1.1 .02 ]' + mvnrnd(0,0.05); % Estimated from document
e_B = e_B/norm(e_B);

e_b_NoiseMag = 0;
e_b_linearNoise = wgn(3,1,e_b_NoiseMag,'linear')';
e_b1 = [ 1 e_b_linearNoise(1) e_b_linearNoise(1) ]';
e_b2 = [ e_b_linearNoise(2) 1 e_b_linearNoise(2) ]';
e_b3 = [ e_b_linearNoise(3) e_b_linearNoise(3) 1 ]';

% Convert modulation direction vectors to unit vectors
e_b1 = e_b1/norm(e_b1)';
e_b2 = e_b2/norm(e_b2)'; 
e_b3 = e_b3/norm(e_b3)'; 

% Create time array 
time = (1/fs):(1/fs):duration;

% External B field guassian white noise 
B_ext_noise = wgn(length(time),1,0,'linear')';

% Modulation pulse frequencies [Hz]
fm1 = 9; 
fm2 = 16;
fm3 = 20;

% Modulation linear guassian white noise [pT]
fm1_noise = wgn(length(time),1,0,'linear')';
fm2_noise = wgn(length(time),1,0,'linear')';
fm3_noise = wgn(length(time),1,0,'linear')';

B_ext = mean(norm(B_external));

% modulation amplitude (ideal case: b1 = b2 = b3 = b)
b1 = 0.5*B_ext*10^(-3);
b2 = B_ext*10^(-3);
b3 = 2*B_ext*10^(-3);

% modulation linear guassian white noise [pT]
b1_noise = wgn(length(time),1,0,'linear')';
b2_noise = wgn(length(time),1,0,'linear')';
b3_noise = wgn(length(time),1,0,'linear')';

% Produce the resulting B field vector components [ x; y; z ]
B_ext_vec = zeros(3,length(time));
b1_vec = zeros(3,length(time));
b2_vec = zeros(3,length(time));
b3_vec = zeros(3,length(time));
for i = 1:3
    % external contribution
    B_ext_vec(i,:) = (B_ext + B_ext_noise)*e_B(i);
    % coil contributions
    b1_vec(i,:) = (b1 + b1_noise).*cos(2*pi.*(fm1 + fm1_noise).*time)*e_b1(i);
    b2_vec(i,:) = (b2 + b2_noise).*cos(2*pi.*(fm2 + fm2_noise).*time)*e_b2(i);
    b3_vec(i,:) = (b3 + b3_noise).*cos(2*pi.*(fm3 + fm3_noise).*time)*e_b3(i);
end


% total B field
B_tot_vec = B_ext_vec + b1_vec + b2_vec + b3_vec;
B_tot = vecnorm(B_tot_vec);
%plot(timeStamp,B_tot)
B_ext_vec = B_ext_vec';

% Spectral analysis of total B field (FFT)
B_tot_fft = fft(B_tot);
B_tot_P2 = abs(B_tot_fft)/(length(B_tot ));
B_tot_P1 = B_tot_P2(1:(length(B_tot )/2) + 1);
B_tot_P1(2:end-1) = 2*B_tot_P1(2:end-1);
B_tot_f = fs*(0:(length(B_tot)/2))/length(B_tot);

end