%% Protocol

%%%
%% General
%%%

%Original pointing first cross scan
N_scan = 5; % number of scan
N_scan_cs = N_scan; %cs : cross scan
N_scan_total = N_scan + N_scan_cs;

%The beta step between two scan (fixed by the protocol Large map, and
%our convention)
scan_step = -348/2; % in arcsec

%number of time sampling point per scan
time_point = 250; % use also (mainly ?) for simulation
N_time_point = ones(N_scan_total,1)*time_point;

%This is the angle of the sensor in the convention defined
%here.
sensor_angle = -42.4*pi/180; %in rad. Certainly only usefull for
                             %simulation. With real data, absolute,
                             %pointing is provided

%This is the angle between the two sensor in each cross scan
cross_scan_angle = 84.8*pi/180; % in rad (fixed by protocol large map)
mat_rot_cross_scan = [cos(cross_scan_angle) -sin(cross_scan_angle); ...
                      sin(cross_scan_angle) cos(cross_scan_angle)];

%%%
%% The speed vector
%%%
speed_norm = 30; % in arcsec*second^-1. Or 60 (fixed by the protocol
                 % Large Map)
speed_angle = 42.4*pi/180; % in rad. Angle between xaxis OF THE SENSOR
                           % and the speed direction. Fixed by
                           % protocol LargeMap
speed_sign = repmat([1 -1], 1, ceil(N_scan/2));
speed_scan = speed_norm*speed_sign(1:N_scan); %(1:N_scan) to not take
                                              %an extra scan if N_scan
                                              %is impair.

%This is zero by our convention. See comment at the top
speed_alpha_scan = speed_scan*cos(0);
speed_beta_scan = speed_scan*sin(0); %in arcsec*s^-1 the_speed is a 2*N_scan
                                     %vector with alpha norm in the first
                                     %line and beta norm in the second line.

%The speed of the first cross scan
the_speeds = [speed_alpha_scan; speed_beta_scan];
%complete with the cross scan (simply apply a rotation)
the_speeds = [the_speeds mat_rot_cross_scan*the_speeds];

%Reshape to have a 2x1xN_scan matrix. First line is alpha, second is
%beta.
the_speeds = reshape(the_speeds,2,1,size(the_speeds,2));

%The four speed vector
unique_speed = [speed_norm , -speed_norm ; 0 , 0 ];
unique_speed = [unique_speed mat_rot_cross_scan*unique_speed];

Nspeed = length(unique_speed(1,:));

