function [H, b, Beta] = sim_data_carolina_withIGRF(B,time)

%time vector
t = length(time);

%% Create helmholtz data
% Three orders of magnitude smaller than external field 

beta_1 = mean(B(:,1))*10^(-2);
beta_2 = mean(B(:,2))*10^(-2);
beta_3 = mean(B(:,3))*10^(-2);

Beta = [beta_1, 0, 0; 0, beta_2,0; 0, 0 beta_3];

%Frequency being driven in each Helmholtz
f_helm_1 = 9; %Hz
f_helm_2 = 16;
f_helm_3 = 25;

%Direction of coils - ideally 1 
%Would need to add noise to this, and eventually some sort of offset 
e_helm_1 = [1 0 0];
e_helm_2 = [0 1 0];
e_helm_3 = [0 0 1];

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

e_helm = eye(3);

B_helm_total1 = [B_helm_1; zeros(1,length(B_helm_1)); zeros(1,length(B_helm_1))];
B_helm_total2 = [zeros(1,length(B_helm_1)); B_helm_2; zeros(1,length(B_helm_1))];
B_helm_total3 = [zeros(1,length(B_helm_1)); zeros(1,length(B_helm_1)); B_helm_3];

%External field 
B_ext_vec = B'+B_helm_total1+B_helm_total2+B_helm_total3;
B_true = B';

%Testing h 
e_1 = repmat([1, 0, 0],length(B),1);
e_2 = repmat([0, 1, 0],length(B),1);
e_3 = repmat([0, 0, 1],length(B),1);

h1 = beta_1*dot(B_true,e_1')/vecnorm(B_true);
h2 = beta_2*dot(B_true,e_2')/vecnorm(B_true);
h3 = beta_3*dot(B_true,e_3')/vecnorm(B_true);

Betas = [beta_1,0,0;0,beta_2,0;0,0,beta_3];
ttest = [h1,h2,h3]*inv(Betas)*mean(vecnorm(B_true));
test_section1 = [h1,h2,h3]*inv(Betas);
test = test_section1'*vecnorm(B_true);

H = [h1,h2,h3];

%Take scalar value 
B_ext_mag = vecnorm(B_true);
    
% pull out DC component
b = mean(vecnorm(B_true));

    

end