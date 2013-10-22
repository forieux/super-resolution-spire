clear all

%% Base
placemount = '/mnt/space/results/'
expname = 'aa2';
system(['mkdir -p ',placemount,expname]);
addpath('../../')
addpath('../../utils')
addpath('../../simulator')
addpath('../../libspire')

randn('state',0)
rand('state',0)

%%
paramsInstrument
paramsObservation; simulatePointing;
paramsSky
physical_sigma_coef = sigma_coef;

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

sigmaInst = 1e4;
true_sigma_coef = physical_sigma_coef - sigmaInst/2;

[Hrond250 Hdirect250] = computeRI(10, 1, band_250, ...
                                  central_wavelength_250, ...
                                  true_sigma_coef, sigma_alpha, ...
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

%% Options
cgoptions.thresold = 1e-8;
cgoptions.maxIter = 40;
%% cgoptions.numfig = 1000;
criterion = 1e-4;

init = zeros(Nalpha, Nbeta);

hypersInit = zeros(2,3);
hypersInit(1,1) = gammaB250;
hypersInit(2,1) = 3.5e+11;

meanInst = physical_sigma_coef;
mean_sigma_coef = meanInst;
riParams = {10, 1, band_250, central_wavelength_250, mean_sigma_coef, ...
            sigma_alpha, sigma_beta, alpha_step, beta_step, Nalpha, ...
            Nbeta, Nspeed, Norder, unique_speed, opt_efficiency, ...
            time_constante, gain};
position = 5;

obs250 = calcObsPix(coefs250, Nalpha, Nbeta, N_scan_total, 20);
%obs250 = ones(size(obs250))

burnin = 200;
maxIter = 2000;
[skyEap gnChain gxChain instChain rate] = myopicUsmse(init, hypersInit, data250, ...
                                       Hrond250, index250, coefs250, ...
                                       offsets, regOp, Nalpha, ...
                                       Nbeta, Norder, N_scan_total, ...
                                       Nbolo250, Nspeed, ...
                                       unique_speed, the_speeds, ...
                                       sigmaInst, meanInst, ...
                                       riParams, position, ...
                                       criterion, burnin, maxIter, ...
                                       cgoptions, 100, obs250);

disp(['Vrai :', num2str(true_sigma_coef)])
disp(['Mean (sigma) :', num2str(meanInst), ' (', num2str(sigmaInst),')'])
disp(['Mean + sigma :', num2str(meanInst + sigmaInst)])
disp(['Est :', num2str(mean(instChain(500:end)))])

name = [placemount,expname,'/doux_mean_above_true_s1e4'];
try
    save(name)
end

