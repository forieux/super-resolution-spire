clear all

%% Base
placemount = '/mnt/space/results/'
expname = 'aa2-CirrusNS'
system(['mkdir -p ',placemount,expname]);
addpath('../../')
addpath('../../utils')
addpath('../../simulator')
addpath('../../libspire')

randn('state',0)
rand('state',0)

%%
paramsInstrument
paramsObservation;
simulatePointing;
paramsSky

%% Precomputation
bound = 4*sigma_coef*max(band_520) + 60; % arcsec
[alpha beta] = computeAxis(pointing250, sky_alpha_period, ...
                           sky_beta_period, ...
                           bound, N_scan_total);

%% Precomputation of index and redondancy
[alphaCoadd, betaCoadd, pointing250coadd, index250coadd, ...
 coefs250coadd, Nalphacoadd, Nbetacoadd] = computeIndex(alpha, beta, ...
                                                        pointing250, ...
                                                        6, 6, ...
                                                        Nbolo250, ...
                                                        N_scan_total);

[alpha, beta, pointing250, index250, coefs250, Nalpha, Nbeta] = ...
    computeIndex(alpha, beta, pointing250, alpha_step, beta_step, ...
                 Nbolo250, N_scan_total);

%% The coaddition without offsets and estimation of offsets
offsets = zeros(1,Nbolo250);

%% Inversion
init = zeros(Nalpha, Nbeta);
objectMean = zeros(size(init));

[Hrond250 Hdirect250] = computeRI(10, 1, band_250, ...
                                  central_wavelength_250, ...
                                  sigma_coef, sigma_alpha, ...
                                  sigma_beta, alpha_step, beta_step, ...
                                  Nalpha, Nbeta, Nspeed, Norder, ...
                                  unique_speed, opt_efficiency, ...
                                  time_constante, gain);
%% Gain
G = Hrond250(1);

[diffAlpha diffBeta circDalpha circDbeta] = computeReg(sigma_alpha, ...
                                                       sigma_beta, ...
                                                       sky_alpha_period, ...
                                                       sky_beta_period, ...
                                                       5, Nalpha, ...
                                                       Nbeta);

regOp = circDalpha + circDbeta;
regOp(1) = 0;

%% Simulate data
std250 = (1e-3); gammaB250 = 1/std250^2;
std360 = std250; gammaB360 = 1/std360^2;
std520 = std360; gammaB520 = 1/std520^2;

std = [std250 std360 std520];
gammaB = [gammaB250, gammaB360, gammaB520];

[data250 output sky] = ...
    simulateData(['../../data&true/simul_ism_SPS_comp'], Nalpha, ...
                 Nbeta, Norder, std250, 10^(-4), Hrond250, ...
                 index250, Nbolo250, Nspeed, N_scan_total, unique_speed, the_speeds);

coaddPSW = dirtymap(data250, index250coadd, zeros(Nbolo250, 1), Nalphacoadd, Nbetacoadd, N_scan_total, Nbolo250);
coadd_interp = imresize(coaddPSW, size(sky), 'nearest');
mesure_sky = sky;
mesure_sky(find(coadd_interp == 0)) = 0;
disp(num2str(sum(abs(coadd_interp(:)/Hrond250(1) - mesure_sky(:)).^2) / sum(abs(mesure_sky(:)).^2)))

%% Options
cgoptions.thresold = 1e-8;
cgoptions.maxIter = 40;

%% cgoptions.numfig = 1000;
the_mean = 0;
for iscan = 1:N_scan_total
   the_mean = the_mean + mean(mean(data250{iscan}));
end

criterion = 1e-4;
burnin = 400; %% If you don't know this value it is more than 100
maxIter = 2500; %% If you don't know it is more 200

init = ones(Nalpha, Nbeta) * the_mean / Hrond250(1,1,1,1);

hypersInit = zeros(2,3);
hypersInit(1,1) = gammaB250;
hypersInit(2,1) = 1e6;

[skyEap gnChain gxChain] = usmse(init, hypersInit, data250, Hrond250, index250, coefs250, offsets, regOp, Nalpha, Nbeta, Norder, N_scan_total, Nbolo250, Nspeed, unique_speed, the_speeds, criterion, burnin, maxIter, cgoptions, 'pla', [placemount,expname]);

mesure_skyEap = skyEap;
mesure_skyEap(find(coadd_interp == 0)) = 0;
disp(num2str(sum(abs(mesure_skyEap(:) - mesure_sky(:)).^2) / sum(abs(mesure_sky(:)).^2)))

save
