% PLOTBEARTH Plot the Earth's magnetic field lines using the IGRF.
% 
% Plots a globe and a number of magnetic field lines starting at each point
% specified by the vectors in lat_start and lon_start. Both distance and
% nsteps should be the same length as lat_start. The plot will spin if spin
% is true and will continue to spin until the user hits CTRL+C.
clear;
close all;
font = 'Times New Roman';
axis_font = 12;
title_font = 12;

addpath('./EKF_code','./igrf_code','./orbit_code');

% *Get the Satellite TLE*
ID = 25544; % ISS
[TLE] = getSatelliteTLE(ID);

% *Convert to Orbital Elements*
[OE] = TLE2OrbitalElements(TLE);
fprintf(1,['Kepler Elements for satelliteID %d epoch %s:\n' ...
    '\ta [km] = %f\n\te = %f\n\ti [deg] = %f\n\tomega [deg] = %f\n' ...
    '\tOmega [deg] = %f\n\tM [deg] = %f\n'], floor(OE.satelliteID), ...
    datestr(OE.epoch),OE.a_km, OE.e, OE.i_deg, OE.omega_deg, ...
    OE.Omega_deg, OE.M_deg);
a_m = OE.a_km*1e3 + 2000*1e3;
e = OE.e;
M_deg = OE.M_deg;

time = datenum([2019 7 17 6 30 0]);
lat_start = 30:15:60; % Geodetic latitudes in degrees.
lon_start = 0:30:330; % Geodetic longitudes in degrees.
alt_start = 0; % Altitude in km.
distance = -sign(lat_start).*[30e3 70e3 150e3]; % km.
nsteps = abs(distance)/10;
spin = false;
% Get the magnetic field line points.
lat = zeros(max(nsteps(:))+1, numel(lat_start)*numel(lon_start));
lon = zeros(max(nsteps(:))+1, numel(lat_start)*numel(lon_start));
alt = zeros(max(nsteps(:))+1, numel(lat_start)*numel(lon_start));
for index1 = 1:numel(lat_start)
    for index2 = 1:numel(lon_start)
        [lat(1:nsteps(index1)+1, ...
            index1*(numel(lon_start)-1)+index2), lon(1:nsteps(index1)+1, ...
            index1*(numel(lon_start)-1)+index2), alt(1:nsteps(index1)+1, ...
            index1*(numel(lon_start)-1)+index2)] = ...
            igrfline(time, lat_start(index1), lon_start(index2), ...
            alt_start, 'geod', distance(index1), nsteps(index1));
    end
end
% Plot the magnetic field lines.
figure;
hold on;
% If the mapping toolbox is not available, use non-mapping toolbox
% functions to plot the globe and lines.
if ~license('test', 'MAP_Toolbox')
    
    % WGS84 parameters.
    a = 6378.137; f = 1/298.257223563;
    b = a*(1 - f); e2 = 1 - (b/a)^2; ep2 = (a/b)^2 - 1;
    
    % Plot a globe.
    load('topo.mat', 'topo', 'topomap1');
    [xe, ye, ze] = ellipsoid(0, 0, 0, a, a, b, 100);
    surface(-xe, -ye, ze, 'FaceColor', 'texture', ...
        'EdgeColor', 'none', 'CData', topo);
    colormap(topomap1);
    
    % Convert lla to xyz.
    [x, y, z] = geod2ecef(lat, lon, alt*1e3); % geod coord
    x = x/1e3; y = y/1e3; z = z/1e3;          % geod coord
    % [x, y, z] = sph2cart(lon*pi/180, lat*pi/180, alt); % geoc coord
    
    % Make the plots.
    plot3(x, y, z, 'r');
    axis equal;
    
% Otherwise, use mapping toolbox functions to plot globe and lines.
else
    load('topo.mat', 'topo', 'topolegend');
    axesm('globe', 'Geoid', 6371.2)
    meshm(topo, topolegend); demcmap(topo);
    % [x, y, z] = sph2cart(lon*pi/180, lat*pi/180, alt*1e3); % geoc coord
    % [lat, lon, alt] = ecef2geod(x, y, z); alt = alt/1e3;   % geoc coord
    plot3m(lat, lon, alt, 'r'); % geod coord
end
% Set the plot background to black.
set(gcf, 'color', 'k');
axis off;
title(['Magnetic Field Lines at ' datestr(time)], 'FontName', font, ...
    'FontSize', title_font, 'Color', 'w');
% Spin the plot indefinitely.
index = 0;
view(mod(index, 360), 23.5); % Earth's axis tilts by 23.5 degrees
while spin
    view(mod(index, 360), 23.5); % Earth's axis tilts by 23.5 degrees
    drawnow;
    pause(0.1);
    index = index - 5;
end

% *Step4_OrbitStateVector2ECR.m*
% *Purpose*
%  To show how to convert the orbital elements to the inertial p-q orbital
%  plane.  The p-axis is through the center of the orbit to perigee.  The q
%  axis is through the focus (center of Earth) and normal to the p axis.


% *Standard Gravitational Parameter*
% The standard gravitational parameter $$\mu$$ of a celestial body is the
% product of the gravitational constant G and the mass M of the body. 
mu = 3.98618e14; % [m3/s2] Earth's geocentric gravitational constant


% *Orbital Plane Coordinates*
%  p_m - [m] coordinate along axis through center and perigee
%  q_m - [m] coordinate passing through focus and perpendicular to p-axis
%  dpdt_m_per_s = [rad/s] p component velocity
%  dqdt_m_per_s = [rad/s] q component velocity
n_rad_per_s = sqrt(mu/a_m^3);  % [rad/s] mean motion
n_deg_per_s = rad2deg(n_rad_per_s); % [deg/s] mean motion
M_rad = deg2rad(M_deg);
E_rad = M_rad; 
dE = 99999;
eps = 1e-6; % [rad] control precision of Newton's method solution
while (abs(dE) > eps)
    dE = (E_rad - e * sin(E_rad) - M_rad)/(1 - e * cos(E_rad));
    E_rad = E_rad -  dE;
end
p_m = a_m*(cos(E_rad) - e);
q_m = a_m*sqrt(1 - e^2)*sin(E_rad);

dMdt_rad_per_s = n_rad_per_s;
dEdt_rad_per_s = dMdt_rad_per_s/(1 - e*cos(E_rad));
dpdt_m_per_s = -a_m*sin(E_rad)*dEdt_rad_per_s;
dqdt_m_per_s = a_m*cos(E_rad)*dEdt_rad_per_s*sqrt(1 - e^2);
E_deg_epoch = rad2deg(E_rad); 

% *Rotate To ECI*
Rz_Omega = [ ...
    [cosd(OE.Omega_deg) sind(OE.Omega_deg) 0]; ...
    [-sind(OE.Omega_deg) cosd(OE.Omega_deg) 0]; ...
    [0 0 1]];
Rx_i = [ ...
    [1 0 0]; ...
    [0 cosd(OE.i_deg) sind(OE.i_deg)]; ...
    [0 -sind(OE.i_deg) cosd(OE.i_deg)]];
Rz_omega = [ ...
    [cosd(OE.omega_deg) sind(OE.omega_deg) 0]; ...
    [-sind(OE.omega_deg) cosd(OE.omega_deg) 0]; ...
    [0 0 1]];

% time of epoch
[Year,Month,Day,H,M,S] = datevec(OE.epoch);
HourUTC = H + M/60.0 + S/3600.0;
jd = juliandate(Year,Month,Day,HourUTC,0,0);
jd0 = juliandate(Year,Month,Day,0,0,0);
% form time in Julian centuries from J2000
T = (jd - 2451545.0d0)./36525.0d0;
% D0 = (jd0 - 2451545.0d0);
% % [deg] GMST = Sidereal Time at Greenwich
% GMST = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*HourUTC + 0.000026*T^2;
% % [deg] Sidereal Time
% ST_deg = 100.46061837 + 36000.770053608*T + 0.000387933*T^2 - T^3/38710000;
% GST_deg = ST_deg + 1.00273790935*HourUTC*15;
% %LST_deg = GST_deg - ObsLongitude_degW;

Theta_deg = 100.460618375 + 36000.770053608336*T + 0.0003879333*T^2 + 15*H + M/4 + mod(S/240,360);

Rz_hour = [ ...
    [cosd(Theta_deg) sind(Theta_deg) 0]; ...
    [-sind(Theta_deg) cosd(Theta_deg) 0]; ...
    [0 0 1]];

% position of satellite at epoch in the orbit pq plane
r_pq = [p_m q_m 0]';
% position of satellite at epoch in ECI coordinates

omega_deg = OE.omega_deg;
Omega_deg = OE.Omega_deg;
i_deg = OE.i_deg;

r_ECI = inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pq;
r_LLA = eci2lla(r_ECI',datevec(datenum(Year, Month, Day, H, M, S)),'IAU-2000/2006');


Evals = 0:1:360.0; % [deg] values of the eccentric anomaly around orbit 
Orbit_p = a_m*(cosd(Evals)-e); % [m] orbit positions
Orbit_q = a_m*sqrt(1 - e^2)*sind(Evals); % [m] orbit positions
deltaT_s = ((Evals-E_deg_epoch) - e*sind(Evals-E_deg_epoch))/n_deg_per_s; % [s] time since epoch along orbit

% only due one at a time... <sigh>
Orbit_ECI = zeros(numel(deltaT_s),3);
Orbit_LLA = zeros(numel(deltaT_s),3);
for ipt = 1:size(Orbit_ECI,1)
    r_pq = [Orbit_p(ipt) Orbit_q(ipt) 0]';
    Orbit_ECI(ipt,:) = [inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pq]'; %[Rz_Omega*Rx_i*Rz_omega*r_pq]';
    lla = eci2lla(Orbit_ECI(ipt,:),datevec(datenum(Year, Month, Day, H, M, S+deltaT_s(ipt))),'IAU-2000/2006');
    Orbit_LLA(ipt,:) = lla;
end

% *Plot Cartesian Coordinates*
hold on;
plot3(Orbit_ECI(:,1)/1000,Orbit_ECI(:,2)/1000,Orbit_ECI(:,3)/1000, 'Linewidth',2.5)
xlabel('ECI x [m]');
ylabel('ECI y [m]');
zlabel('ECI z [m]');
title('Satellite Orbit in ECI Coordinates');
grid on

% *Plot Map*
figure('color','white');
earth = imread('ear0xuu2.jpg');
lv= size(earth,1);
lh= size(earth,2);
lats =  (1:lv)*180/lv - 90;
lons =  (1:lh)*360/lh - 180;
image(lons, -lats, earth)
hold on;
set(gca,'ydir','normal');
grid on
plot(Orbit_LLA(:,2),Orbit_LLA(:,1),'.r');
plot(r_LLA(2),r_LLA(1),'p','MarkerFaceColor',[1 0 0.7],'MarkerEdgeColor','none','Markersize',14);
set(gca,'XTick',[-180:30:180]);
set(gca,'YTick',[-90:30:90]);
set(gca,'Xcolor',0.3*ones(1,3));
set(gca,'Ycolor',0.3*ones(1,3));
title(['Ground Track for Epoch ' datestr(OE.epoch)])
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
axis equal
axis tight
set(gca,'fontweight','bold');


[latitude, longitude, altitude] = ecef2geod(Orbit_ECI(:,1), Orbit_ECI(:,2), Orbit_ECI(:,3), 1e-3);

% Call to IGRF for a series of locations
[Bx, By, Bz] = igrf(time, latitude, longitude, altitude, 'geocentric');

nSteps = length(latitude);

% Add perturbations according to noiseModel and pertMagnitude
%noiseModel    = 'gaussian';
%noiseModel    = 'students-t';
%noiseModel    = 'gmm';
noiseModel    = 'exp';
%noiseModel    = 'none';

pertMagnitude = 1e-5;

if strcmp(noiseModel, 'gaussian')
    % Gaussian noise
    pert = randn(3,nSteps)*pertMagnitude;
elseif strcmp(noiseModel, 'students-t')
    % Student's T distribution noise
    pert = trnd(3,nSteps)*pertMagnitude;
elseif strcmp(noiseModel, 'gmm')
    % Gaussian mixture model noise
    mu = randn(5, 1)*10;
    gm = gmdistribution(mu,diag(pertMagnitude));
    pert = reshape(gm.random(3*nSteps)*pertMagnitude, 3, nSteps);
elseif strcmp(noiseModel, 'none')
    % No noise added
    pert =  zeros(3,nSteps);
elseif strcmp(noiseModel, 'exp')
    % Two-sided exponential noise
    pert = -sign(randn(3,nSteps))*pertMagnitude.*log(rand(3,nSteps));
end

Bx = Bx + pert(1,:)';
By = By + pert(2,:)';
Bz = Bz + 10*pert(3,:)';

% Plot B-field components
figure(); subplot(3,1,1); plot(Bx); ylabel('B_x [nT]'); grid on;
subplot(3,1,2); plot(By); ylabel('B_y [nT]'); grid on;
subplot(3,1,3); plot(Bz); ylabel('B_z [nT]'); grid on; xlabel('Sample');
