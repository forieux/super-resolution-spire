%% This script contains the physical parameter of the instrument.

%% Optics

mirorDiameter = 3.285; % meter
hollDiameter = 0.56; % meter
focal = 8.68; % meter

opt_efficiency = 0.95*0.9; % wihtout dim, between 0 and 1

%% Filter

%% Taken from sensitivity model Griffin etal.
band_250 = [213, 290]*10^-6; % in meter
band_360 = [301, 405]*10^-6; % in meter
band_520 = [410, 611]*10^-6; % in meter
                                                             
%% This is the real central wavelength (with our above definition of the
%% filter bound)
central_wavelength_250 = mean(band_250,2); % in meter
central_wavelength_360 = mean(band_360,2); % in meter
central_wavelength_520 = mean(band_520,2); % in meter

%% Feedhorns 

%% Gaussian width
sigma_coef = 24.6/(353*10^-6*sqrt(log(256)));
%% arcsec*meter^-1 (because the wavelength is in meter., and I want the
%% width in arcsec..); For the value see Ferlet07. The origine of
%% sqrt(log256) is that the FHWM of a gaussian is sqrt(log256)*sigma

%% Bolometer. Linear model

time_constante = 19.45e-3; % s 
gain = 3.31e8; % V/W

%% Lecture electronique
temporal_sampling_periode = 1/15; % in s

%% New time constante
time_constante = time_constante + 0; % s for 5Hz filter

%% Number of bolometer
Nbolo250 = size(data{1}{1},2);
Nbolo360 = size(data{2}{1},2);
Nbolo520 = size(data{3}{1},2);

%% Corelated noise parameters (empirical, fixed after lecture of bruce
%% simulator manual, see slide etc)

frequencyCut = 1e-3; % Hz in absolute frequency
frequencyCutE = frequencyCut*temporal_sampling_periode; % reduced frequency
inclination = 2;
stdC = 6e-1;
