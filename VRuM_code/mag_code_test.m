
clear; clf;

nDataPoints = 1000;
pertMagnitude = 1e-3;

% Noise models to add to IGRF
%noiseModel    = 'gaussian';   % AWGN
%noiseModel    = 'gmm';        % Gaussian mixture model 
%noiseModel    = 'exp';        % exponential noise
noiseModel    = 'students-t'; % student's-t noise
%noiseModel    = 'none';       % no noise

[B_total, B_external] = generateVRuMdataFromIGRF(nDataPoints, ...
                                                 noiseModel, ...
                                                 pertMagnitude);
                                                                                                  
figure(1); subplot(2,1,1); hold on; grid on;
plot(B_external', '--')
plot(B_total', '.');
legend('B_{x,ext}','B_{y,ext}','B_{z,ext}',...
       'B_x','B_y','B_z');
xlabel('Data Points')
ylabel('Magnetic Field Intensity [nT]')

subplot(2,1,2); hold on; grid on;
histogram(B_total(1,:) - B_external(1,:))
histogram(B_total(2,:) - B_external(2,:))
histogram(B_total(3,:) - B_external(3,:))
