function [H,B_ext_mag, b, B] = sim_data_carolina_withIGRF(noiseModel,B,time,f)

%time vector
t = length(time);

%% Create helmholtz data
% Three orders of magnitude smaller than external field 

beta_1 = 50;
beta_2 = 25;
beta_3 = 70;

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

B_helm_total1 = [B_helm_1; zeros(1,length(B_helm_1)); zeros(1,length(B_helm_1))];
B_helm_total2 = [zeros(1,length(B_helm_1)); B_helm_2; zeros(1,length(B_helm_1))];
B_helm_total3 = [zeros(1,length(B_helm_1)); zeros(1,length(B_helm_1)); B_helm_3];

%External field 
B_ext_vec = B+B_helm_total1+B_helm_total2+B_helm_total3;

B_ext_mag = vecnorm(B_ext_vec);


N = length(vecnorm(B)); 
xdft = fft(vecnorm(B));
xdft = xdft(1:(N-1)/2+1);
psdx = (1/(N))*abs(xdft);
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:100/N:100/2;

figure;
semilogy(freq,psdx);

[H_peak,fm] = maxk(psdx(2:end),3);

[f_sort,ind] = sort(fm);
H = H_peak(ind);
    
% pull out DC component
b = psdx(1,1);

figure;
semilogy(freq,psdx);
set(gca,'YScale','log','YMinorTick','on')
    title('Spectra of B Field Magnitude')
    

end