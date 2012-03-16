clear all

%% Base
% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

%%

placemount = '/espace/orieux/results/'
expname = 'variousDirty'
system(['mkdir -p ',placemount,expname]);

%% 
paramsInstrument
paramsObservation; simulatePointing;
paramsSky

%% Precomputation
[alpha, beta, start_end_position, pointing250, pointing360, pointing520, ...
 index250, index360, index520, coefs250, coefs360, coefs520] = ...
    computeIndex(alpha, beta, start_end_position, pointing250, pointing360, ...
                 pointing520, alpha_step, beta_step);

[alpha_mm, beta_mm, start_end_position_mm, pointing250_mm, pointing360_mm, ...
 pointing520_mm, index250_mm, index360_mm, index520_mm, coefs250_mm, ...
 coefs360_mm, coefs520_mm] = computeIndex(alpha, beta, start_end_position, ...
                                          pointing250, pointing360, ...
                                          pointing520, 6, 6);

precalculsInvariantRI

%% Data simulation
std250 = (1e-3); gammaB250 = 1/std250^2;
std360 = std250; gammaB360 = 1/std360^2;
std520 = std360; gammaB520 = 1/std520^2;

std = [std250 std360 std520];
gammaB = [gammaB250, gammaB360, gammaB520];

bands = [1 1 1];

[data output data250 data360 data520 sky] = simulateData('../data&true/simul_ism_SPS_comp', ...
                                              Nalpha, Nbeta, Norder, [1 1 ...
                    1], std, 10^(-4), Hrond250, Hrond360, Hrond520, index250, ...
                                              index360, index520);

%% Retropojection of data and TF of it (dirty map)
map = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                   coefs360, coefs520, Nalpha, Nbeta);

madmap = dirtymap(data, bands, index250_mm, index360_mm, index520_mm, ...
                  coefs250_mm, coefs360_mm, coefs520_mm, length(unique(alpha_mm)) ...
                  , length(unique(beta_mm)));

mapmean = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                       coefs360, coefs520, Nalpha, Nbeta, 'm');

mapnan = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                  coefs360, coefs520, Nalpha, Nbeta, 'n');

save([placemount,expname,'/dirtymaps'])

