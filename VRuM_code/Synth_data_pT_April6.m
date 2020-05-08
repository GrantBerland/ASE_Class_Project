function [B_tot_P1, B_ext_vec] = Synth_data_pT_April6(temp_noise,plot_flag)
%% User Input Parameters

% Set sampling freq., duration
fs = 100; % Hz ( 100 Hz in Swarm doc )
duration = 0.3; % sec

% EXTERNAL B FIELD
% guassian white noise in external B field (linear [pT])
B_ext_linearNoise = 0 + temp_noise; 
% B_ext_linearNoise = 1500 + temp_noise; % pT (Noise at Fs of 100 Hz from previous EMI test: 1500 nt)

% external B field direction vector
e_B_noiseMag = 0;
e_B_linearNoise = mvnrnd(0,0.05);
e_B = [ .8  1.1 .02 ]' + e_B_linearNoise; % Estimated from document

% MODULATION B FIELDS
% guassian white noise in modulation pulsation (linear [Hz])
% noise here seems to have the largest effect
%{
fm1_linearNoise = (1 + temp_noise/10)*10^-12; % Hz
fm2_linearNoise = (1 + temp_noise/10)*10^-12; % Hz
fm3_linearNoise = (1 + temp_noise/10)*10^-12; % Hz
%}
fm1_linearNoise = 0; % Hz
fm2_linearNoise = 0; % Hz
fm3_linearNoise = 0; % Hz

% guassian white noise in modulation amplitude (linear [pT])
%{
b1_linearNoise = 200 + temp_noise/10; % pT
b2_linearNoise = 200 + temp_noise/10; % pT
b3_linearNoise = 200 + temp_noise/10; % pT
%}
b1_linearNoise = 15000; % pT
b2_linearNoise = 15000; % pT
b3_linearNoise = 15000; % pT

% modulation direction vectors (ideal case: aligned with coord. sys.)
e_b_NoiseMag = 0;
e_b_linearNoise = wgn(3,1,e_b_NoiseMag,'linear')';
e_b1 = [ 1 e_b_linearNoise(1) e_b_linearNoise(1) ]';
e_b2 = [ e_b_linearNoise(2) 1 e_b_linearNoise(2) ]';
e_b3 = [ e_b_linearNoise(3) e_b_linearNoise(3) 1 ]';

%% Synthetic Data Production

% create time stamp 
timeStamp = (1/fs):(1/fs):duration;

% external B field amplitude [pT]
B_ext = 50000*10^3;
% external B field guassian white noise 
B_ext_noise = wgn(length(timeStamp),1,B_ext_linearNoise,'linear')'; 
% convert external B field direction vector to unit vector
temp = norm(e_B);
e_B = [ e_B(1)/temp e_B(2)/temp e_B(3)/temp ]'; 

% modulation pulsations [Hz]
fm1 = 9; 
fm2 = 16;
fm3 = 20;
% modulation linear guassian white noise [pT]
fm1_noise = wgn(length(timeStamp),1,fm1_linearNoise,'linear')';
fm2_noise = wgn(length(timeStamp),1,fm2_linearNoise,'linear')';
fm3_noise = wgn(length(timeStamp),1,fm3_linearNoise,'linear')';

% modulation amplitude (ideal case: b1 = b2 = b3 = b)
b1 = 0.5*B_ext*10^(-3);
b2 = B_ext*10^(-3);
b3 = 2*B_ext*10^(-3);
% modulation linear guassian white noise [pT]
b1_noise = wgn(length(timeStamp),1,b1_linearNoise,'linear')';
b2_noise = wgn(length(timeStamp),1,b2_linearNoise,'linear')';
b3_noise = wgn(length(timeStamp),1,b3_linearNoise,'linear')';
% convert modulation direction vectors to unit vectors
temp = norm(e_b1);
e_b1 = [ e_b1(1)/temp e_b1(2)/temp e_b1(3)/temp ]';
temp = norm(e_b2);
e_b2 = [ e_b2(1)/temp e_b2(2)/temp e_b2(3)/temp ]'; 
temp = norm(e_b3);
e_b3 = [ e_b3(1)/temp e_b3(2)/temp e_b3(3)/temp ]'; 

% Produce the resulting B field vector components [ x; y; z ]
B_ext_vec = zeros(3,length(timeStamp));
b1_vec = zeros(3,length(timeStamp));
b2_vec = zeros(3,length(timeStamp));
b3_vec = zeros(3,length(timeStamp));
for i = 1:3
    % external contribution
    B_ext_vec(i,:) = (B_ext + B_ext_noise)*e_B(i);
    % coil contributions
    b1_vec(i,1:length(timeStamp)) = (b1 + b1_noise).*cos(2*pi.*(fm1 + fm1_noise).*timeStamp)*e_b1(i);
    b2_vec(i,1:length(timeStamp)) = (b2 + b2_noise).*cos(2*pi.*(fm2 + fm2_noise).*timeStamp)*e_b2(i);
    b3_vec(i,1:length(timeStamp)) = (b3 + b3_noise).*cos(2*pi.*(fm3 + fm3_noise).*timeStamp)*e_b3(i);
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

%% Plotting

if plot_flag
    plot(B_tot_f,B_tot_P1);
    set(gca,'YScale','log','YMinorTick','on');
    title('Spectra of B Field Magnitude')
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [pT/\surd Hz]')
    xlim([-1,B_tot_f(end)])
    ylim([.8, inf])
end

end