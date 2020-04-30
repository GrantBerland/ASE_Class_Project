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

%% Retrive synthetic data

fs = 100;
time = 1/fs:1/fs:1;

time_B = [0,15, 30, 45,60];
%B_points = [50000,50020,50050,50080,50110,50095,50105,50100,50080,50040,49990,49960,50010,49975,50010,50050];
B_points = [50000 49980 50010 49985 49996];
time_want = 0:1/fs:60;
B = interp1(time_B,B_points,time_want,'spline');

j = 1;
noise = 0;
b_mag = zeros(1,59);

for i = 1:59
    [H, b, Ext_vec ,Beta] = sim_data_carolina_perdata_set(B(j:j+100),noise,time_want(1:101),fs);
 %Save true field 
True_ext(:,j:j+100) = Ext_vec;
b_mag(i) = b;

j = j+101;
%% trying to figure out what SVD is suppossed to be 
lam = Beta;
A_true = eye(3);
G_true_inv = lam*A_true*lam';
G_true = inv(G_true_inv);


%% Determine B_curly (projection of B_ext_vec onto e_b directions) and B_ext_vec

B_cur = b.*H*inv(lam);
B_ext(:,i) = B_cur*inv(A_true);
  

end
%%
B_ext_points = B_ext;
B_ext_times = 0.5:1:59;
B_ext_x = interp1(B_ext_times,B_ext_points(1,:),time_want,'spline');
B_ext_y = interp1(B_ext_times,B_ext_points(2,:),time_want,'spline');
B_ext_z = interp1(B_ext_times,B_ext_points(3,:),time_want,'spline');

h = figure;
set(h,'Color',[1 1 1]);
set(h,'DefaultAxesFontSize',15);
set(h,'DefaultAxesFontName', 'Times New Roman')
set(h,'DefaultTextFontname', 'Times New Roman')
plot(time_want,B_ext_x,'LineWidth',2)
hold on
plot(time_want(:,1:end-42),True_ext(1,:),'LineWidth',2)
title('B_x true magnetic field vs VRuM measured field')
legend('VRuM Field x component','True Field x component')
ylabel('Magnetic Field (nT)')
xlabel('time')
grid on
ylim([31420,31880])

h = figure;
set(h,'Color',[1 1 1]);
set(h,'DefaultAxesFontSize',15);
set(h,'DefaultAxesFontName', 'Times New Roman')
set(h,'DefaultTextFontname', 'Times New Roman')
plot(time_want,B_ext_y,'LineWidth',2)
hold on
plot(time_want(:,1:end-42),True_ext(2,:),'LineWidth',2)
title('B_y true magnetic field vs VRuM measured field')
legend('VRuM Field y component','True Field y component')
ylabel('Magnetic Field (nT)')
xlabel('time')
grid on
ylim([33420,33880])

h = figure;
set(h,'Color',[1 1 1]);
set(h,'DefaultAxesFontSize',15);
set(h,'DefaultAxesFontName', 'Times New Roman')
set(h,'DefaultTextFontname', 'Times New Roman')
plot(time_want,B_ext_z,'LineWidth',2)
hold on
plot(time_want(:,1:end-42),True_ext(3,:),'LineWidth',2)
title('B_z true magnetic field vs VRuM measured field')
legend('VRuM Field z component','True Field z component')
ylabel('Magnetic Field (nT)')
xlabel('time')
grid on
ylim([19680,19970])
