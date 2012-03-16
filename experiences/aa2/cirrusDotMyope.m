clear all

%% Base
placemount = '/mnt/space/results/'
expname = 'cirrusDotMyope'
system(['mkdir -p ',placemount,expname]);
addpath(genpath('../../'))

randn('state',0)
rand('state',0)

%%
paramsInstrument
paramsObservation; simulatePointing;
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

%% Simulate data
std250 = (1e-3); gammaB250 = 1/std250^2;
std360 = std250; gammaB360 = 1/std360^2;
std520 = std360; gammaB520 = 1/std520^2;

std = [std250 std360 std520];
gammaB = [gammaB250, gammaB360, gammaB520];

[data250 output sky] = ...
    simulateData(['../../data&true/simul_cond_ism_SPS_comp'], Nalpha, ...
                 Nbeta, Norder, std250, 10^(-4), Hrond250, ...
                 index250, Nbolo250, Nspeed, N_scan_total, unique_speed, the_speeds);

%% Options
cgoptions.thresold = 1e-8;
cgoptions.maxIter = 40;
%% cgoptions.numfig = 1000;
criterion = 1e-4;
burnin = 400; %% If you don't know this value it is more than 100
maxIter = 3000; %% If you don't know it is more 200

init = zeros(Nalpha, Nbeta);

hypersInit = zeros(2,3);
hypersInit(1,1) = gammaB250;
hypersInit(2,1) = 1e6;

riParams = {10, 1, band_250, central_wavelength_250, sigma_coef, ...
            sigma_alpha, sigma_beta, alpha_step, beta_step, Nalpha, ...
            Nbeta, Nspeed, Norder, unique_speed, opt_efficiency, ...
            time_constante, gain};

sigmaInst = 1e3;
sigma_coef
meanInst = sigma_coef + sigmaInst
position = 5;

tic
[skyEap gnChain gxChain instChain rate] = myopicUsmse(init, hypersInit, data250, ...
                                       Hrond250, index250, coefs250, ...
                                       offsets, regOp, Nalpha, ...
                                       Nbeta, Norder, N_scan_total, ...
                                       Nbolo250, Nspeed, ...
                                       unique_speed, the_speeds, ...
                                       sigmaInst, meanInst, ...
                                       riParams, position, ...
                                       criterion, burnin, maxIter, ...
                                       cgoptions);

toc

name = [placemount,expname,'/eap'];
try
    save(name)
end


