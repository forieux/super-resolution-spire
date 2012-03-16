%% Protocol

N_scan = 5;
N_scan_total = 10;
N_time_point = [672 671 666 672 666 672 671 672 672 672];

%% This is the angle of the sensor in the convention defined here.
%% Certainly only usefull for simulation. With real data, absolute,
%% pointing is provided
sensor_angle = -42.4*pi/180; %in rad.

%% This is the angle between the two sensor in each cross scan
cross_scan_angle = 84.8*pi/180; % in rad (fixed by protocol large map)
mat_rot_cross_scan = [cos(cross_scan_angle) -sin(cross_scan_angle); ...
                      sin(cross_scan_angle) cos(cross_scan_angle)];

%% The speed vector
speed_norm = 30; % 30 or 60 in arcsec*second^-1
%% Angle between xaxis OF THE SENSOR and the speed direction. Fixed by
%% protocol LargeMap
speed_angle = 42.4*pi/180; % in rad.
speed_sign = repmat([1 -1], 1, ceil(N_scan/2));
%% To not take an extra scan if N_scan is impair.
speed_scan = speed_norm*speed_sign(1:N_scan); %(1:N_scan)

%% This is zero by our convention. See comment at the top
speed_alpha_scan = speed_scan*cos(0);
speed_beta_scan = speed_scan*sin(0); % arcsec*s^-1 

%% the_speed is a 2*N_scan vector with alpha norm in the first line and
%% beta norm in the second line. The speed of the first cross scan
the_speeds = [speed_alpha_scan; speed_beta_scan];
%% complete with the cross scan (simply apply a rotation)
the_speeds = [the_speeds mat_rot_cross_scan*the_speeds];

%% Reshape to have a 2x1xN_scan matrix. First line is alpha, second is
%% beta.
the_speeds = reshape(the_speeds,2,1,size(the_speeds,2));

unique_speed = [speed_norm , -speed_norm ; 0 , 0 ];
unique_speed = [unique_speed mat_rot_cross_scan*unique_speed];

Nspeed = length(unique_speed(1,:));

disp('The speed vector is not good !')
