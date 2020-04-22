function [B_ext_mag, Ext_vec] = sim_data_carolina(temp_noise,time,f)

%time vector
t = length(time);

%% Creating external field with noise
%Will end up getting this for IGRF for adv state project 
%External Field nano tesla
Ext_field = 50000; 

%Make it for 3 vectors
B_ext_noise_vec = temp_noise;

%For now will have the same noise for all three axis
%Divide the magnitude of 50 000 nT into components
%For now will make the magnitude divide pretty evenly

direction_B_ext = [0.63; 0.67; 0.395];

Ext_vec_noiseless = Ext_field*direction_B_ext/norm(direction_B_ext);
Ext_vec_noiseless_time = Ext_vec_noiseless.*ones(3,t);
Ext_vec = Ext_vec_noiseless_time+B_ext_noise_vec; 

%% Create helmholtz data
% Three orders of magnitude smaller than external field 

beta_1 = Ext_vec_noiseless(1)*10^(-3);
beta_2 = Ext_vec_noiseless(2)*10^(-3);
beta_3 = Ext_vec_noiseless(3)*10^(-3);

%Frequency being driven in each Helmholtz
f_helm_1 = 9; %Hz
f_helm_2 = 16;
f_helm_3 = 25;

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
B_helm_1 = beta_1*cos(2*pi*f_helm_1*time);
B_helm_2 = beta_2*cos(2*pi*f_helm_2*time);
B_helm_3 = beta_3*cos(2*pi*f_helm_3*time);

%Direction vectors of helmholtz coils, for now they are perfectly
%ortogonal
e_helm = [1,0,0;0,1,0;0,0,1];
e_helm = e_helm/norm(e_helm);

B_helm_total = [B_helm_1',B_helm_2',B_helm_3']*e_helm;


%External field 
B_ext_vec(1,:) = Ext_vec(1,:) + B_helm_1;
B_ext_vec(2,:) = Ext_vec(2,:) + B_helm_2;
B_ext_vec(3,:) = Ext_vec(3,:) + B_helm_3;




%Magnitude of external field - ie what magnetometer sees

B_ext_mag = vecnorm(B_ext_vec);

N = length(B_ext_mag);

xdft = fft(B_ext_mag);
xdft = xdft(1:N/2+1);
psdx = (1/(N))*abs(xdft);
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:f/N:f/2;

figure;
semilogy(freq,psdx);
grid on





end