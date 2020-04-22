function [B_ext_mag, B] = sim_data_carolina_withIGRF(noiseModel,time,f)

%time vector
t = length(time);

%% Creating external field with noise
%Coming from IGRF model
% Noise models to add to IGRF

B = generate_B_field_dynamics(noiseModel,1);

%% Create helmholtz data
% Three orders of magnitude smaller than external field 

beta_1 = 50;
beta_2 = 50;
beta_3 = 50;

%Frequency being driven in each Helmholtz
f_helm_1 = 9; %Hz
f_helm_2 = 16;
f_helm_3 = 24;

%Direction of coils - ideally 1 
%Would need to add noise to this, and eventually some sort of offset 
%Do these have to vecnorm to 1?
e_helm_1 = 1;
e_helm_2 = 1;
e_helm_3 = 1;

%Initialize vector of Helmholtz
B_helm_1 = zeros(length(time),1);
B_helm_2 = zeros(length(time),1);
B_helm_3 = zeros(length(time),1);
%Magnetic Field made by each helmhotz - for now noiseless
B_helm_1 = beta_1*cos(2*pi*f_helm_1*time)+randn(1);
B_helm_2 = beta_2*cos(2*pi*f_helm_2*time)+randn(1);
B_helm_3 = beta_3*cos(2*pi*f_helm_3*time)+randn(1);

%External field 
B_ext_x = abs(B(:,1)) + B_helm_1';
B_ext_y = abs(B(:,2)) + B_helm_2';
B_ext_z = abs(B(:,3)) + B_helm_3';

%Magnitude of external field - ie what magnetometer sees

B_ext_vec = [B_ext_x,B_ext_y,B_ext_z];

B_ext_mag = vecnorm(B_ext_vec');

B_tot_fft = fft(B_ext_mag);
B_tot_P2 = abs(B_tot_fft)/1000;
B_tot_P1 = B_tot_P2(1:(1001-1)/2 + 1);
B_tot_P1(2:end-1) = 2*B_tot_P1(2:end-1);
B_tot_f = f*(0:1000/2)/(1001);

figure;
semilogy(B_tot_f,B_tot_P1);
set(gca,'YScale','log','YMinorTick','on')
    title('Spectra of B Field Magnitude')


end