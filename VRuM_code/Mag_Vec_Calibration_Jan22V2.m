%% House Keeping
clear all
close all
clc

addpath('../igrf_code')
%% Set Up
% enter number of data sets
N = 1000;
% pre-allocate variables for speed
H = zeros(N,3);
fm = zeros(N,3);
b = zeros(N,1);
S_inv = zeros(3,3);
B_curly_vec = zeros(N,3);
B_ext_vec = zeros(N,3);
rng(100)
%% Retrive synthetic data

[data, B_ext_truel] = generateVRuMdataFromIGRF(N, 'none', 5e-4);

% pull out principle harmonics ( find peaks may not be the best way)
[H(i,:), fm(i,:)] = maxk(data(i,2:end),3); % endpoints excluded

b(i,1) = data(i,1);

%% Retirive G matrix entries
% Implament SVD algorithm
[U,S,V] = svd(H,'econ');
s = svd(H,'econ');
% if condition > 10^12 (with doubles) than A is ill-conditioned
condition = max(diag(S))/min(diag(S));

% set the threshold for the SVD 
t = 1e-10;
% set small singular values (<t) to zero
S_inv = inv(S);
S_inv(S_inv < t) = 0;


result = eye(length(H));


% Solve the system for G
G = V*S_inv*U'*result*U*S_inv*V';
% Take the inverse of G_mat usign built in matlab function
G_inv = inv(G);


%test
AA = H*G*H';
%% Determine alpha, theta, gamma, beta(1,2, and 3) See Appindix A in Swarm

% define diagonal lam values
lam_11 = sqrt(G_inv(1,1)); % beta1
lam_22 = sqrt(G_inv(2,2)); % beta2
lam_33 = sqrt(G_inv(3,3)); % beta3  We need to be careful because this isn't guaranteed to be positive

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
   B_curly_vec(i,:) =  b(i,1).* H(i,:) / (Lambda_mat);
   B_ext_vec(i,:) =  B_curly_vec(i,:) / (A_mat);
end

% Difference between true and calculated external B field
B_true_mean = mean(B_ext_true);
B_ext_mean = mean(B_ext_vec);
B_error = B_ext_true - B_ext_vec;
%%

time = linspace(0,10,1000); h = figure;
set(h,'Color',[1 1 1]);
set(h,'DefaultAxesFontSize',15);
set(h,'DefaultAxesFontName', 'Times New Roman')
set(h,'DefaultTextFontname', 'Times New Roman')

subplot(3,1,1); hold on;
plot(time',B_ext_true(:,1)*0.001,'LineWidth',2);
plot(time',B_ext_vec(:,1)*0.001,'LineWidth',2);
title('Components of True Magnetic Field vs Estimated Magnetic field over time')
ylabel('x - B field (nT)'); legend('True B Field','Estimated B Field')

subplot(3,1,2); hold on;
plot(time',B_ext_true(:,2)*0.001,'LineWidth',2);
plot(time',B_ext_vec(:,2)*0.001,'LineWidth',2);
ylabel('y - B field (nT)'); legend('True B Field','Estimated B Field')

subplot(3,1,3); hold on;
plot(time',B_ext_true(:,3)*0.001,'LineWidth',2);
plot(time',B_ext_vec(:,3)*0.001,'LineWidth',2);
xlabel('Time (seconds)')
ylabel('z - B field (nT)'); legend('True B Field','Estimated B Field')
