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
data_sets = 40;
fs = 100;
time_end = 600;
time = 1/fs:1/fs:time_end;

%preallocate matrices
data = zeros(1,length(time));


%% Retrive synthetic data

pertMagnitude = 0.5;

%Works right now - for eq 2 - 8
%B_extfield_overtime = [50000,35000,10000,-5000,5000,-8000,20000,32000,48000];
B_extfield_overtime = [70000,45000,10000,-5000,15000,-8000,20000,32000,50000];
%Testing my data
j = 1;
j_add = 99;
for i = 1:9
    
    Ext_field = B_extfield_overtime(i);
    noiseModel = randn(3,1+j_add)*pertMagnitude;
    [H(i,:),b(i),B_true] = sim_data_carolina_perdata_set(Ext_field,noiseModel,time(j:j+j_add),fs);
    j = j+j_add;
end

close all

%% Retirive G matrix entries

%% trying to figure out what SVD is suppossed to be 
lam = [31.4713,0,0;0,33.4695,0;0,0,19.7320];
A_true = eye(3);
G_true_inv = lam*A_true*lam';
G_true = inv(G_true_inv);

%does this equal I
shouldbeI = H*G_true*H';

%Solving the equations individually 
[G_indv] = Mag_Vec_solveindividually(H);

G_ind_inv = inv(G_indv);


%% Using SVD
% Implament SVD algorithm
[U,S,V] = svd(H,'econ');

%The matrix on the right side of the system of linear equations HGH = 1 for
%k runs - see Gravrand et al pg 954
result = eye(3);

shouldbeII = ones(3,3);
% Solve the system for G - see attached SVD document 
G = V*inv(S)*U'*result*U*inv(S)*V';

% Take the inverse of G_mat usign built in matlab function
G_inv = inv(G);

%Up to this point something is wrong!
%test - AA should equal identity matrix 
AA = H*G*H';



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


for i = 1:9
   B_curly_vec(i,:) =  b(1,i).* H(i,:)* inv(Lambda_mat);
   B_ext_vec(i,:) =  B_curly_vec(i,:)*inv(A_mat);
end
