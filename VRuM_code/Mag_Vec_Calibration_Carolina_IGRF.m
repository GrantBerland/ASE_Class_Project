%% Magnetic Field vectorization Method - March version
% Script that find the vector magnetic field of the external field the 
% magnetometer is trying to measure by analyzing the spectra of the 
% magnetic field magnitude 
% Based on Gravrand et al., 2001, "On the calibration of a vectorial 4He
% pumped magnetometer"

%Carolina Pena and Jaykob Velazquez

clear all
close all
%clc
rng(100)

%addpath('./EKF_code','./igrf_code','./orbit_code');

%% Set Up

% enter number data sets
data_sets = 5;
fs = 100;
time_end = 20;
time = 1/fs:1/fs:20;

%preallocate matrices
data = zeros(1,length(time));
B_true = zeros(length(time),3);


%% Retrive synthetic data

pertMagnitude = 0.5;

noiseModel = randn(3,length(time))*pertMagnitude;

%Testing my data
[data,B_true] = sim_data_carolina_withIGRF(noiseModel,time,fs);
%%
%Find spectogram of data for window of 0.4 seconds
h = figure; 
    set(h,'Color',[1 1 1]);
    set(h,'DefaultAxesFontSize',30)
[s,t,f] = spectrogram(data,100,[],[],fs);
colormap(jet)
close all
s = abs(s);

freq = mean(s');



%%
%Generate fft
j = 1;
j_add = length(time)/data_sets-1;
for i = 1:data_sets
    [fft_data(i,:),B_f] = generate_fft(data(j:j+j_add),fs);
    j = j+j_add;
    
     % pull out principle harmonics and the frequencies of these peaks 
    [H(i,:),fm(i,:)] = maxk(fft_data(i,2:end),3); % endpoints excluded
    
    %The H (ie magnitude of peaks) that are found above are in the order of
    %largest peak to lowest peak, but that is not the order of the
    %frequecies - so need to sort the H vector such that it follows the
    %order of the freuquencies
    [f_sort(i,:),ind(i,:)] = sort(fm(i,:));
    H(i,:) = H(i,ind(i,:));
    
    % pull out DC component
    b(i,1) = fft_data(i,1);
end

%% Retirive G matrix entries
% Implament SVD algorithm
[U,S,V] = svd(H,'econ');

%The matrix on the right side of the system of linear equations HGH = 1 for
%k runs - see Gravrand et al pg 954
result = eye(length(H(:,1)));

% Solve the system for G - see attached SVD document 
G = V*inv(S)*U'*result*U*inv(S)*V';

% Take the inverse of G_mat usign built in matlab function
G_inv = inv(G);

%Up to this point something is wrong!
%test - AA should equal identity matrix 
AA = H(1:3,1:3)*G*H(1:3,1:3)';


%% Determine alpha, theta, gamma, beta(1,2, and 3) See Appindix A in Gravrand et al

% define diagonal lam values
lam_11 = sqrt(G_inv(1,1)); % beta1
lam_22 = sqrt(G_inv(2,2)); % beta2
lam_33 = sqrt(G_inv(3,3)); % beta3  


% determine alpha
sin_alpha = -G_inv(1,2)/(lam_11*lam_22);
alpha = asin(sin_alpha);
cos_alpha = cos(alpha);
%cos_alpha = sqrt(1-(sin_alpha)^2);

% determine p and q
p = G_inv(1,3)/(lam_11*lam_33);
q = (G_inv(2,3)/(lam_22*lam_33)+(sin_alpha*p))/(cos_alpha);

% determine theta and gamma
tan_theta_sq = p^2/(1-p^2-q^2);
tan_gamma_sq = q^2/(1-p^2-q^2);

%% Define Lambda, C, and A matricies

% Lambda matrix
Lambda_mat = [ lam_11 0 0; 0 lam_22 0; 0 0 lam_33 ]; 

% C matrix
C_33 = 1/(sqrt(1+tan_theta_sq+tan_gamma_sq)); 
C_mat = [ 1 0 0; -sin_alpha cos_alpha 0; p q C_33 ];

% A matrix
A_mat = C_mat*C_mat';

%% Determine B_curly (projection of B_ext_vec onto e_b directions) and B_ext_vec


for i = 1:data_sets
   B_curly_vec(i,:) =  b(i,1).* H(i,:)* inv(Lambda_mat);
   B_ext_vec(i,:) =  B_curly_vec(i,:)*inv(A_mat);
end
