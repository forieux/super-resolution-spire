%	$Id: paramsInstrument.m,v 1.7.1.4 2009/08/31 11:43:02 orieux Exp $	
%
% Time-stamp: <2010-07-06 15:39:06 orieux>

%This script contains the physical parameter of the instrument.

%% Optics

%%% Airy
mirorDiameter = 3.285; % meter
hollDiameter = 0.56; % meter
focal = 8.68; % meter

opt_efficiency = 0.95*0.9; % wihtout dim, between 0 and 1

%% Filter
central_wavelength_250 = 250*10^-6; % in meter
central_wavelength_360 = 360*10^-6; % in meter
central_wavelength_520 = 520*10^-6; % in meter

band_250 = [250 - (360 - 250)/2, 250 + (360 - 250)/2]*10^-6; % in meter
band_360 = [360 - (360 - 250)/2, 360 + (520 - 360)/2]*10^-6; % in meter
band_520 = [520 - (520 - 360)/2, 520 + (520 - 360)/2]*10^-6; % in meter

% Taken from sensitivity model Griffin etal.
band_250 = [213, 290]*10^-6; % in meter
band_360 = [301, 405]*10^-6; % in meter
band_520 = [410, 611]*10^-6; % in meter
                                                             
%This is the real central wavelength (with our above definition of
%the filter bound)
central_wavelength_250 = mean(band_250,2); % in meter
central_wavelength_360 = mean(band_360,2); % in meter
central_wavelength_520 = mean(band_520,2); % in meter

%% Feedhorns 

%%% Gaussian width

%%% Spatial sampling step in the direction alpha (pa) and beta (pb)
p_alpha_250 = 2*central_wavelength_250/mirorDiameter; % radian
                                                      % (small angle hyp)
p_alpha_250 = rad2arcsec(p_alpha_250); % arcsec
p_beta_250  = p_alpha_250/sqrt(3);   % arcsec

p_alpha_360 = 2*central_wavelength_360/mirorDiameter; % radian
                                                      % (small angle hyp)
p_alpha_360 = rad2arcsec(p_alpha_360); % arcsec
p_beta_360  = p_alpha_360/sqrt(3);   % arcsec

p_alpha_520 = 2*central_wavelength_520/mirorDiameter; % radian
                                                      % (small angle hyp)
p_alpha_520 = rad2arcsec(p_alpha_520); % arcsec
p_beta_520  = p_alpha_520/sqrt(3);   % arcsec

%% Complete Gaussian optics
% sigma_coef = 15/(350*10^-6); % arcsec*meter^-1 (because the wavelength is
%                              % in meter., and I want the width in
%                              % arcsec..); For the value see Ferlet07.

%% BE CAREFUL : the true is most this one. The previous is for simulation !
sigma_coef = 24.6/(353*10^-6*sqrt(log(256))); 
%arcsec*meter^-1 (because the wavelength is in meter., and I want the width
%in arcsec..); For the value see Ferlet07. The origine of sqrt(log256) is
%that the FHWM of a gaussian is sqrt(log256)*sigma

%% Bolometer

%%% Thermal model. This parameter are not use actually. Since the
%model is linear, we only need the gain and the time constante.

% polarisation = 0.015; % Volt
% load_resistance = 16000; % Ohm

% gain_zero = 65e-12; % Watt/Kelvin
% temperature_zero = 0.2; % Kelvin
% beta_link = 1.7;

% resistance_star = 80; % Ohm
% temperature_star = 41; % Kelvin
% n_resistance = 0.5; % this parameter is without dimention

%%% Linear model

time_constante = 19.45e-3; % s 
gain = 3.31e8; % V/W

% time_constante = 0.75; % for simulation
% gain = 1/time_constante; % for simulation to have a normalized response

%% Lecture electronique
temporal_sampling_periode = 1/15; % in s

time_constante = time_constante + 0.2; % s 


%% Number of bolometer
Nbolo250 = 139;
Nbolo360 = 88;
Nbolo520 = 43;

%% Corelated noise parameters (empirical, fixed after lecture of bruce
%simulator manual, see slide etc)

frequencyCut = 1e-3; % Hz in absolute frequency
frequencyCutE = frequencyCut*temporal_sampling_periode; % reduced frequency
inclination = 2;
stdC = 6e-1;


%	$Log: paramsInstrument.m,v $
%	Revision 1.7.1.4  2009/08/31 11:43:02  orieux
%	Correct width
%
%	Revision 1.7.1.3  2009/08/31 11:37:54  orieux
%	Correct width
%
%	Revision 1.7.1.2  2009/08/26 08:37:26  orieux
%	Correlated signal paramters add
%
%	Revision 1.7  2009/08/05 14:59:10  orieux
%	*** empty log message ***
%
%	Revision 1.6  2009/06/01 08:24:12  orieux
%	Gain = 1/tau to have a normalized response
%
%	Revision 1.5  2009/05/20 12:17:35  orieux
%	Add the number of bolometer
%
%	Revision 1.4  2009/05/15 14:39:01  orieux
%	*** empty log message ***
%
%	Revision 1.3  2009/05/15 14:02:31  orieux
%	Plus d'affichage à chaque itération
%
%	Revision 1.2  2009/05/15 13:11:11  orieux
%	*** empty log message ***
%	
