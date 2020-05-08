%% Magnetic Field vectorization Method - March version
% Script that find the vector magnetic field of the external field the 
% magnetometer is trying to measure by analyzing the spectra of the 
% magnetic field magnitude 
% Based on Gravrand et al., 2001, "On the calibration of a vectorial 4He
% pumped magnetometer"

%Carolina Pena and Jaykob Velazquez

clear all
close all
clc
rng(100)

%% Set Up
% enter number of data sets
N = 30;
% pre-allocate variables for speed
H = zeros(N,3);
fm = zeros(N,3);
b = zeros(N,1);
S_inv = zeros(3,3);
B_curly_vec = zeros(N,3);
B_ext_vec = zeros(N,3);

%Set rng so the random numbers do not change in each run
rng(100)


%Make noise
for j = 1:N
    random_noise(j) = abs(randn*200);
end 

%Gather data

for i = 1:N
    %Get one of the random noise values
   noise = random_noise(i);
   %Gather data - scalar data and true data 
   [spectral_data(i,:),B_ext_true] = Synth_data_pT_April6(noise,false);  
   
   %Take out 
   b(:,1) = spectral_data(i,1);
   
   %Get the peaks and frequencies
   [H_unordered(i,:),H_index(i,:)] = maxk(spectral_data(i,2:end),3);
   
    %The H (ie magnitude of peaks) that are found above are in the order of
    %largest peak to lowest peak, but that is not the order of the
    %frequecies - so need to sort the H vector such that it follows the
    %order of the freuquencies
    [f_sort(i,:),index(i,:)] = sort(H_index(i,:));
    H(i,:) = H_unordered(i,index(i,:));   
end

%Get SVD

[U,S,V] = svd(H);

%S_inv = inv(S(1:3,1:3));

[U_n,S_n,V_n] = mysvd(H);

S_inv = inv(S_n);

G_n = V_n*S_inv*U_n'*eye(30)*U_n*S_inv*V_n';

Test = H*G_n*H';

G_inv = inv(G_n);

% define diagonal lam values
lam_11 = sqrt(G_inv(1,1)); % beta1
lam_22 = sqrt(G_inv(2,2)); % beta2
lam_33 = sqrt(G_inv(3,3)); % beta3  

%Calculating B field

% determine alpha
sin_alpha = -G_inv(1,2)/(lam_11*lam_22);
cos_alpha = sqrt(1-(sin_alpha^2));

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
for i = 1:N
   B_curly_vec(i,:) =  b(i,1).* H(i,:)* inv(Lambda_mat);
   B_ext_vec(i,:) =  B_curly_vec(i,:)*inv(A_mat);
end
%Difference between true and calculated external B field

%%

time = linspace(0,10,1000);
h = figure;
set(h,'Color',[1 1 1]);
set(h,'DefaultAxesFontSize',15);
set(h,'DefaultAxesFontName', 'Times New Roman')
set(h,'DefaultTextFontname', 'Times New Roman')
subplot(3,1,1)
plot(B_ext_true(:,1)*0.001,'LineWidth',2);
hold on
plot(B_ext_vec(:,1)*0.001,'LineWidth',2);
title('Components of True Magnetic Field vs Estimated Magnetic field over time')
hold off
ylabel('x - B field (nT)')
legend('True B Field','Estimated B Field')
subplot(3,1,2)
plot(B_ext_true(:,2)*0.001,'LineWidth',2);
hold on
plot(B_ext_vec(:,2)*0.001,'LineWidth',2);
ylabel('y - B field (nT)')
legend('True B Field','Estimated B Field')
hold off
subplot(3,1,3)
plot(B_ext_true(:,3)*0.001,'LineWidth',2);
hold on
plot(B_ext_vec(:,3)*0.001,'LineWidth',2);
xlabel('Time (seconds)')
ylabel('z - B field (nT)')
legend('True B Field','Estimated B Field')
hold off

