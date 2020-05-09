function measurements = Mag_Vec_Calibration_Carolina_shortcut_IGRF(noiseModel,permag,plot_status,initial_state,dt,nMeasurements)

% Script that find the vector magnetic field of the external field the 
% magnetometer is trying to measure by analyzing the spectra of the 
% magnetic field magnitude 
% Based on Gravrand et al., 2001, "On the calibration of a vectorial 4He
% pumped magnetometer"


rng(100)


%% Set Up
fs = 100;
time_end = 1;
time_perdata = 1/fs:1/fs:0.05;
total_time = 1/fs:1/fs:1;

b_pertime = zeros(1,total_time(end)-1);

%True B field 

B = generate_B_field_measurement(noiseModel, 0, 0, initial_state(1:3), initial_state(4:6), dt, nMeasurements);

%% Retrive synthetic data

j = 1;
for i = 1:24
    [H, b, Beta] = sim_data_carolina_withIGRF(B(j:j+4,:),time_perdata);
    j = j+4;
%Saving B mag as a vector 
b_pertime(i) = b;

%% trying to figure out what SVD is suppossed to be 
%Save Beta as a diagonal matrix
lam = Beta;
%A_true is indentity matrix because we are assuming coils are perpendicular
A_true = eye(3);

%% Determine B_curly (projection of B_ext_vec onto e_b directions) and B_ext_vec

B_cur = b.*H*inv(lam);
B_ext(:,i) = B_cur*inv(A_true);
end 
time_want = 1/fs:1/fs:1;
B_ext_points = B_ext;
B_ext_times = 0.04:0.04:1-1/fs;
B_ext_x = interp1(B_ext_times,B_ext_points(1,:),time_want,'spline');
B_ext_y = interp1(B_ext_times,B_ext_points(2,:),time_want,'spline');
B_ext_z = interp1(B_ext_times,B_ext_points(3,:),time_want,'spline');

measurements = [B_ext_x',B_ext_y',B_ext_z'];

%{
h = figure;
set(h,'Color',[1 1 1]);
set(h,'DefaultAxesFontSize',15);
set(h,'DefaultAxesFontName', 'Times New Roman')
set(h,'DefaultTextFontname', 'Times New Roman')
plot(time_want,B_ext_x,'*')
hold on
plot(time_want,B(:,1))
title('B_x true magnetic field vs VRuM measured field')
legend('VRuM Field x component','True Field x component')
ylabel('Magnetic Field (nT)')
xlabel('time')
grid on

h = figure;
set(h,'Color',[1 1 1]);
set(h,'DefaultAxesFontSize',15);
set(h,'DefaultAxesFontName', 'Times New Roman')
set(h,'DefaultTextFontname', 'Times New Roman')
plot(time_want,B_ext_y,'*')
hold on
plot(time_want,B(:,2))
title('B_y true magnetic field vs VRuM measured field')
legend('VRuM Field y component','True Field y component')
ylabel('Magnetic Field (nT)')
xlabel('time')
grid on

h = figure;
set(h,'Color',[1 1 1]);
set(h,'DefaultAxesFontSize',15);
set(h,'DefaultAxesFontName', 'Times New Roman')
set(h,'DefaultTextFontname', 'Times New Roman')
plot(time_want,B_ext_z,'*')
hold on
plot(time_want,B(:,3))
title('B_z true magnetic field vs VRuM measured field')
legend('VRuM Field z component','True Field z component')
ylabel('Magnetic Field (nT)')
xlabel('time')
grid on

%}

end