function [H,b, External_true, Beta] = sim_data_carolina_perdata_set(Ext_field,temp_noise,time,f)

%time vector
t = length(time);

%% Creating external field with noise
%Will end up getting this for IGRF for adv state project 
%External Field nano tesla
%Ext_field = 50000; 

%Make it for 3 vectors
B_ext_noise_vec = temp_noise;

%For now will have the same noise for all three axis
%Divide the magnitude of 50 000 nT into components
%For now will make the magnitude divide pretty evenly


direction_B_ext = [0.63; 0.67; 0.395];

Ext_vec = Ext_field.*direction_B_ext/norm(direction_B_ext);
%Ext_vec_noiseless_time = Ext_vec_noiseless.*ones(3,t);
%Ext_vec = Ext_vec_noiseless_time+B_ext_noise_vec; 



%% Create helmholtz data
% Three orders of magnitude smaller than external field 

beta_1 = 6;
beta_2 = 8;
beta_3 = 3;

Beta = [beta_1, 0, 0; 0, beta_2,0; 0, 0 beta_3];

%Frequency being driven in each Helmholtz
f_helm_1 = 9; %Hz
f_helm_2 = 16;
f_helm_3 = 25;

%Direction of coils - ideally 1 
%Would need to add noise to this, and eventually some sort of offset 
%Do these have to vecnorm to 1?
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
%for t = 1:length(B_helm_1)
%    B_helm_total(:,t) = [B_helm_1%(t)',B_helm_2(t)',B_helm_3(t)'] * e_helm;
%end
%B_helm_total = [B_helm_1',B_helm_2',B_helm_3'] * repmat(e_helm, length(B_helm_1), 1);

%B_helm_total1 = B_helm_1' * repmat(e_helm_1,length(B_helm_1),1);
B_helm_total1 = [B_helm_1', zeros(length(B_helm_1), 1),zeros(length(B_helm_1), 1)];
B_helm_total2 = [zeros(length(B_helm_2), 1), B_helm_2', zeros(length(B_helm_2), 1)];
B_helm_total3 = [zeros(length(B_helm_1), 1), zeros(length(B_helm_1), 1), B_helm_3'];

%B_helm_total3 = [0 * B_helm_1; 0.05 * B_helm_2; 1 *  B_helm_3'];


%External field 
Ext_vec = Ext_vec';
B_ext_vec = Ext_vec+B_helm_total1+B_helm_total2+B_helm_total3;

External_true = Ext_vec';


%Magnitude of external field - ie what magnetometer sees


%Testing h 

e_1 = repmat([1, 0, 0],101,1);
e_2 = repmat([0, 1, 0],101,1);
e_3 = repmat([0, 0, 1],101,1);

h1 = beta_1*dot(External_true,e_1')/vecnorm(External_true);
h2 = beta_2*dot(External_true,e_2')/vecnorm(External_true);
h3 = beta_3*dot(External_true,e_3')/vecnorm(External_true);

Betas = [beta_1,0,0;0,beta_2,0;0,0,beta_3];
test = zeros(3,1);
test = [h1,h2,h3]*inv(Betas)*mean(vecnorm(External_true));

B_ext_mag = vecnorm(B_ext_vec');
H = [h1,h2,h3];

N = length(B_ext_mag);

xdft = fft(B_ext_mag);
xdft = xdft(1:(N-1)/2+1);
psdx = (1/(N-1))*abs(xdft);
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:f/(N-1):f/2;

[H_peak,fm] = maxk(psdx(5:end),3);

[f_sort,ind] = sort(fm);
HH = H_peak(ind);
 
%H = [psdx(10), psdx(17), psdx(26)];
% pull out DC component
b = psdx(1,1);

%figure;
%semilogy(freq,psdx);
%grid on




end