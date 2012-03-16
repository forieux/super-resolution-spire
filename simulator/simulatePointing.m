%% Certainly only for simulation. Compute the pointing.

%% The start position of each scan. They are computed here for
%% simulation need. They can be provided

%% start_end_position is a 2xN matrix. The first line is the alpha
%% position of the center of the sensor at the start then the end of
%% each scan. The second line is beta at the start and the end.

start_end_position = [0, N_time_point(1)* temporal_sampling_periode* ...
                    speed_alpha_scan(1); scan_step, scan_step];

for scan = 2:N_scan
  %% compute for the first cross scan. It is simply a shift in beta
  start_end_position_add = [start_end_position(1,(scan-1)*2), ...
                            start_end_position(1,(scan-1)*2) + ...
                            N_time_point(scan)* ...
                            temporal_sampling_periode*speed_alpha_scan(scan); ...
                            scan*scan_step, scan*scan_step];
    
  %% add it to the full vector
  start_end_position = [start_end_position start_end_position_add];
end
%% start_end_position(2,:) = start_end_position(2,:) + scan_step/4;

%% The position for cross scan are computed with the previous position,
%% a rotation of the sensor and a shift_cs. The rotation is normaly
%% fixed by the protocol so for the simulation we apply this angle. It
%% could be different in reallity and it is possible to take into
%% account this info. At the end, what is necessary is the position of
%% the bolometer at each time. How this position is obtained is not
%% importante.
start_end_position_cross_scan = mat_rot_cross_scan* start_end_position;
%% The shift is arbitrary here for simulation.

%% First move to the origine
start_end_position_cross_scan(1,:) = start_end_position_cross_scan(1,:) - ...
    start_end_position_cross_scan(1,1);
start_end_position_cross_scan(2,:) = start_end_position_cross_scan(2,:) - ...
    start_end_position_cross_scan(2,1);

shift_cs = [30 3*scan_step]'; % [alpha beta]' Arbitrary
shift_cs = [start_end_position(1,end-1) start_end_position(2,end-1)]'; % [alpha beta]' Arbitrary
start_end_position_cross_scan = start_end_position_cross_scan + ...
    repmat(shift_cs,1,size(start_end_position,2));

%% Concatenation of the two to get all the start_end position
start_end_position = [start_end_position start_end_position_cross_scan];

%% This a 2x2xN tab. First line is alpha, second is beta. First col
%% is start, second is end. Third dim is the scan.
start_end_position = reshape(start_end_position, [2 2 N_scan_total]);

clear start_end_position_cross_scan start_end_position_add     

%%
%%  The position of all bolometer in each scan
%%

%% Position in sensor system coordinate

%% 250
l_250 = ([0:4]  - 2)*p_alpha_250;
m_250 = ([0:14] - 7)*p_beta_250;

%% The submatrix
l_250_2 = ([0:3]  - 1.5)*p_alpha_250;
m_250_2 = ([0:15] - 7.5)*p_beta_250;

%% 360
l_360 = ([0:3] - 1.5)*p_alpha_360;
m_360 = ([0:12] - 6)*p_beta_360;

l_360_2 = ([0:2] - 1)*p_alpha_360;
m_360_2 = ([0:11] - 5.5)*p_beta_360;

%% 520
l_520 = ([0:2] - 1)*p_alpha_520;
m_520 = ([0:8] - 4)*p_beta_520;

l_520_2 = ([0:1] - 0.5)*p_alpha_520;
m_520_2 = ([0:7] - 3.5)*p_beta_520;

%% Same but Gridded
[PA_250 PB_250] = ndgrid(l_250, m_250);
[PA_250_2 PB_250_2] = ndgrid(l_250_2, m_250_2);
                    
[PA_360 PB_360] = ndgrid(l_360, m_360);
[PA_360_2 PB_360_2] = ndgrid(l_360_2, m_360_2);
                    
[PA_520 PB_520] = ndgrid(l_520, m_520);
[PA_520_2 PB_520_2] = ndgrid(l_520_2, m_520_2);


%% The position in sky coordinate

%% 250

%%  COORD_... is similar to a 2D tab but with two element in each cell: the
%%  alpha and beta coord:
%%  |----------------+----------------|--
%%  | (alpha1,beta1) | (alpha1,beta2) |
%%  |----------------+----------------|--
%%  | (alpha2,beta1) | (alpha2,beta2) |
%%  |----------------+----------------|--
COORD_250 = zeros(2,length(l_250),length(m_250));
COORD_250(1,:,:) = PA_250;
COORD_250(2,:,:) = PB_250;

%%  Each (alpha, beta) is rearanged to a matrix with first line is alpha
%%  and 2 line is beta. We do this to apply easily rotation
%%  |-------------+-------------+---
%%  | alpha_bolo1 | alpha_bolo2 |   
%%  |-------------+-------------+---
%%  | beta_bolo1  | beta_bolo2  |   
%%  |-------------+-------------+---
COORD_250 = reshape(COORD_250,2,prod(size(PA_250)));

%% The submatrix 250_2
COORD_250_2 = zeros(2,length(l_250_2),length(m_250_2));
COORD_250_2(1,:,:) = PA_250_2;
COORD_250_2(2,:,:) = PB_250_2;
COORD_250_2 = reshape(COORD_250_2,2,prod(size(PA_250_2)));
COORD_250 = [COORD_250 COORD_250_2];

%% assigne a unique number to each bolo. This is their name (a number)
%% to identify them.
name_250 = [1:size(COORD_250,2)];

%% 360

COORD_360 = zeros(2,length(l_360),length(m_360));
COORD_360(1,:,:) = PA_360;
COORD_360(2,:,:) = PB_360;
COORD_360 = reshape(COORD_360,2,prod(size(PA_360)));

%% 360_2
COORD_360_2 = zeros(2,length(l_360_2),length(m_360_2));
COORD_360_2(1,:,:) = PA_360_2;
COORD_360_2(2,:,:) = PB_360_2;
COORD_360_2 = reshape(COORD_360_2,2,prod(size(PA_360_2)));
COORD_360 = [COORD_360 COORD_360_2];

%%  assigne a unique number to each bolo
name_360 = [1:size(COORD_360,2)];

%% 520

COORD_520 = zeros(2,length(l_520),length(m_520));
COORD_520(1,:,:) = PA_520;
COORD_520(2,:,:) = PB_520;
COORD_520 = reshape(COORD_520,2,prod(size(PA_520)));

%% 520_2
COORD_520_2 = zeros(2,length(l_520_2),length(m_520_2));
COORD_520_2(1,:,:) = PA_520_2;
COORD_520_2(2,:,:) = PB_520_2;
COORD_520_2 = reshape(COORD_520_2,2,prod(size(PA_520_2)));
COORD_520 = [COORD_520 COORD_520_2];

%%  assigne a unique number to each bolo
name_520 = [1:size(COORD_520,2)];

%% Now we have to compute in coordinate defined here, with the two
%% submatrix concatened.

%% The rotation matrix to apply.
sensor_mat_rot = [cos(sensor_angle) -sin(sensor_angle); sin(sensor_angle) ...
                  cos(sensor_angle)];

%% Repeat in a third dim for the scan (the initial position of each
%% bolo can be different,  a minima because there is crossscan)
COORD_250_scan1 = repmat(sensor_mat_rot*[COORD_250 COORD_250_2],[1 1 N_scan]);
COORD_360_scan1 = repmat(sensor_mat_rot*[COORD_360 COORD_360_2],[1 1 N_scan]);
COORD_520_scan1 = repmat(sensor_mat_rot*[COORD_520 COORD_520_2],[1 1 N_scan]);

%% for the cross scan (the angle is differente)
sensor_mat_rot_cs = [cos(-sensor_angle) -sin(-sensor_angle); ...
                    sin(-sensor_angle) cos(-sensor_angle)];
%%  !! No more angle between cs !!
%% sensor_mat_rot_cs = [1 0; 0 1];

COORD_250_scan2 = repmat(sensor_mat_rot_cs*[COORD_250 COORD_250_2],[1 1 N_scan_cs]);
COORD_360_scan2 = repmat(sensor_mat_rot_cs*[COORD_360 COORD_360_2],[1 1 N_scan_cs]);
COORD_520_scan2 = repmat(sensor_mat_rot_cs*[COORD_520 COORD_520_2],[1 1 N_scan_cs]);

%% allocate before concatenation
COORD_250_scan = zeros(size(COORD_250_scan1,1),size(COORD_250_scan1,2),N_scan_total);
COORD_360_scan = zeros(size(COORD_360_scan1,1),size(COORD_360_scan1,2),N_scan_total);
COORD_520_scan = zeros(size(COORD_520_scan1,1),size(COORD_520_scan1,2),N_scan_total);

%% concatenation; X_scan is 2x1xN_scan tab. First line is alpha
%% position of each bolo without dep. Second line is beta. Thirs dim
%% is the scan.
COORD_250_scan(:,:,1:N_scan) = COORD_250_scan1;
COORD_250_scan(:,:,N_scan+1:N_scan_total) = COORD_250_scan2;
COORD_360_scan(:,:,1:N_scan) = COORD_360_scan1;
COORD_360_scan(:,:,N_scan+1:N_scan_total) = COORD_360_scan2;
COORD_520_scan(:,:,1:N_scan) = COORD_520_scan1;
COORD_520_scan(:,:,N_scan+1:N_scan_total) = COORD_520_scan2;

clear l_250 l_520_2 m_250 m_250_2 l_360 l_360_2 m_360 m_360_2 l_520 ...
    l_520_2 m_520 m_520_2
clear PA_250 PB_250 PA_250_2 PB_250_2 PA_360 PB_360 PA_360_2 PB_360_2 ...
    PA_520 PB_520 PA_520_2 PB_520_2 
clear COORD_250_2 COORD_360_2 COORD_520_2 
clear COORD_250_scan2 COORD_360_scan2 COORD_520_scan2 COORD_250_scan1 ...
    COORD_360_scan1 COORD_520_scan1

%% Positin during scanning

%% Allocation for each point: multiplication by the number of time
%% sample point. 1e ligne is alpha, 2e is beta. Number of column is the
%% number of bolometer. Third dim is scan. We do this to apply easily
%% rotation

for iscan = 1:N_scan_total
    
  %% COORD_X_total are cells (use {}) and not tab. This is because the
  %% number of point in each scan can be (in theory) different. Each
  %% cell is a 2x(Nbr_bolometerxNbr_point). The number of cells is
  %% N_scan_total.
    
  COORD_250_tmp = zeros(2, size(COORD_250,2)*N_time_point(iscan));
  COORD_360_tmp = zeros(2, size(COORD_360,2)*N_time_point(iscan));
  COORD_520_tmp = zeros(2, size(COORD_520,2)*N_time_point(iscan));
    
  dep_alpha = [0:N_time_point(iscan)-1]*the_speeds(1,1,iscan)* ...
      temporal_sampling_periode;
    
  dep_beta = [0:N_time_point(iscan)-1]*the_speeds(2,1,iscan)* ...
      temporal_sampling_periode;
    
  dep = [dep_alpha ; dep_beta];

  %% For each bolo, add all the position
  for n = 0:size(COORD_250,2) - 1
    COORD_250_tmp(1,1+n*N_time_point(iscan):(n+1)* N_time_point(iscan)) = start_end_position(1,1,iscan) + COORD_250_scan(1,n+1,iscan) + dep_alpha;
    COORD_250_tmp(2,1+n*N_time_point(iscan):(n+1)* N_time_point(iscan)) = start_end_position(2,1,iscan) + COORD_250_scan(2,n+1,iscan) + dep_beta;
  end
  pointing250(iscan) = {COORD_250_tmp};
    
  for n = 0:size(COORD_360,2) - 1
    COORD_360_tmp(1,1+n*N_time_point(iscan):(n+1)* N_time_point(iscan)) ...
        = start_end_position(1,1,iscan) + COORD_360_scan(1,n+1,iscan) ...
        + dep_alpha;
    COORD_360_tmp(2,1+n*N_time_point(iscan):(n+1)* N_time_point(iscan)) ...
        = start_end_position(2,1,iscan) + COORD_360_scan(2,n+1,iscan) ...
        + dep_beta;
  end
  pointing360(iscan) = {COORD_360_tmp};
  
  for n = 0:size(COORD_520,2) - 1
    COORD_520_tmp(1,1+n*N_time_point(iscan):(n+1)* N_time_point(iscan)) ...
        = start_end_position(1,1,iscan) + COORD_520_scan(1,n+1,iscan) ...
        + dep_alpha;
    COORD_520_tmp(2,1+n*N_time_point(iscan):(n+1)* N_time_point(iscan)) ...
        = start_end_position(2,1,iscan) + COORD_520_scan(2,n+1,iscan) ...
        + dep_beta;
  end
  pointing520(iscan) = {COORD_520_tmp};
  
end

%%  Each COORD_..._total are now cells with

%%  |----------------+----------------+-----+----------------+--
%%  | bolo1_n1_alpha | bolo1_n2_alpha | ... | bolo2_n1_alpha |  
%%  |----------------+----------------+-----+----------------+--
%%  | bolo1_n1_beta  | bolo1_n2_beta  | ... | bolo2_n2_beta  |  
%%  |----------------+----------------+-----+----------------+--

%%  in each cell. The number of cell is the number of scan

