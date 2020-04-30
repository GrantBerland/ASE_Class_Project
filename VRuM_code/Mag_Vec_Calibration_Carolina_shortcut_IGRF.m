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


%% Set Up

% enter number data sets
data_sets = 40;
fs = 100;
nSteps = 100;
time_end = 600;
time = 1/fs:1/fs:time_end;

%preallocate matrices
data = zeros(1,length(time));

%Real data from IGFR

%Get true data from dynamics
B = generate_B_field_dynamics('gaussian',1);
B = B';

%% Retrive synthetic data

fs = 100;
time = 1/fs:1/fs:1;

j = 1;
for i = 1:3
    [H, b, Beta] = sim_data_carolina_withIGRF('gaussian',B(:,j:j+32),time(j:j+32),fs);
    j = j+33;
end 

%% trying to figure out what SVD is suppossed to be 
lam = Beta;
A_true = eye(3);
G_true_inv = lam*A_true*lam';
G_true = inv(G_true_inv);


%% Determine B_curly (projection of B_ext_vec onto e_b directions) and B_ext_vec

B_cur = b.*H*inv(lam);
B_ext = B_cur*inv(A_true);

for i = 1:9
   B_curly_vec(i,:) =  b(1,i).* H(i,:)* inv(lam);
   B_ext(i,:) =  B_curly_vec(i,:)*inv(A_trye);
end
