clear all

format('long');
%% set path, specially gpac
path(path,'../')
path(path,'../libspire')
path(path,'../utils')
path(path,'../')

racine = '../spire/';
expname = 'ngc/';
system(['mkdir -p ', racine, expname]);
name = '/pv';

load ../data&true/data_fits_image
load ../data&true/pos_fits_image

%% Observation
paramsInstrumentNGC7023_me

data250 = data{1};
N_scan_total = size(data{1},2);
Nbolo250 = size(data{1}{1},2);

pointing250 = calcProj(pointing{1}, N_scan_total);
[unique_speed, the_speeds, Nspeed] = calcSpeed(pointing250, N_scan_total, Nbolo250, temporal_sampling_periode);

%% Sky init
alpha_step = 6; beta_step = 6; Norder = 1;
sky_alpha_period = alpha_step;
sky_beta_period = sky_alpha_period;
sigma_alpha = 0.6*sky_alpha_period;
sigma_beta  = 0.6*sky_beta_period;

bound = 4*sigma_coef*max(band_520) + 4*60; % arcsec
[alpha beta] = computeAxis(pointing250, sky_alpha_period, sky_beta_period, bound, N_scan_total);

%% Precomputation of index and redondancy
[alphaCoadd, betaCoadd, pointing250coadd, index250Coadd, coefs250Coadd, NalphaCoadd, NbetaCoadd] = computeIndex(alpha, beta, pointing250, 6, 6, Nbolo250, N_scan_total);
[alpha, beta, pointing250inv, index250, coefs250, Nalpha, Nbeta] = computeIndex(alpha, beta, pointing250, alpha_step, beta_step, Nbolo250, N_scan_total);

%% Inversion
init = zeros(Nalpha, Nbeta);
objectMean = zeros(size(init));

%% Instrument model
[Hrond250 Hdirect250] = computeRI(10, 1, band_250, central_wavelength_250, sigma_coef, sigma_alpha, sigma_beta, alpha_step, beta_step, Nalpha, Nbeta, Nspeed, Norder, unique_speed, opt_efficiency, time_constante, gain);

break
G = Hrond250(1); % Gain

%% Regularity measure
[diffAlpha diffBeta circDalpha circDbeta] = computeReg(sigma_alpha, sigma_beta, sky_alpha_period, sky_beta_period, 5, Nalpha, Nbeta);
regOp = circDalpha + circDbeta;

%% Conjugate gradient options
cgoptions.thresold = 1e-10; cgoptions.maxIter = 40; %cgoptions.numfig = 1000;

hyp_tab = 5e-7;
offsets = zeros(1,Nbolo250);

keyboard

for hyp_ind = 1:size(hyp_tab,2)
    
    disp(['Hyper ',num2str(hyp_ind),'/',num2str(size(hyp_tab,2))])

    [mapCg offsetsCg ooptimCg] = inversionOffsets(init, cgoptions, data250, [hyp_tab(hyp_ind); 1e1], Hrond250, index250, coefs250, offsets, regOp, objectMean, Nalpha, Nbeta, Norder, N_scan_total, Nbolo250, Nspeed, unique_speed, the_speeds, 5);

end

%% The coaddition without offsets and estimation of offsets
[coaddPSW offsetsPSW] = dirtymapOffsets(data250, offsetsCg, index250Coadd, NalphaCoadd, NbetaCoadd, N_scan_total, Nbolo250, 1); % 20

%%
figure(5)
theAlpha = -100;
ligne = find(alpha <= theAlpha, 1, 'last' );
ligneCoadd = find(alphaCoadd <= theAlpha, 1, 'last' );
clf
plot(beta, mapCg(ligne,:)*G,'r')
hold on
plot(betaCoadd, coaddPSW(ligneCoadd,:))
