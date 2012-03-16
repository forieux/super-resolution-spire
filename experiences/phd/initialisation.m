%% VARIANCECIRRUSDCT - Compute posterior variance of cirrus 

placemount = '/espace/orieux/results/'
expname = '';
system(['mkdir -p ',placemount,expname]);

%% Base
% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

%% 
paramsInstrument
paramsObservation; simulatePointing;
paramsSky

Norder = 2;

%% Precomputation
precalculsInvariantPointing
precalculsInvariantRI

%% Data simulation
std250 = (1e-3); gammaB250 = 1/std250^2;
std360 = std250; gammaB360 = 1/std360^2;
std520 = std360; gammaB520 = 1/std520^2;

std = [std250 std360 std520];
gammaB = [gammaB250, gammaB360, gammaB520];

[data output data250 data360 data520] = simulateData(['../data&true/' ...
                    'simul_ism_SPS_comp'], Nalpha, Nbeta, Norder, [1 1 1], ...
                                                  std, 10^(-4), Hrond250, ...
                                                  Hrond360, Hrond520, ...
                                                  index250, index360, ...
                                                  index520);

name = [placemount,expname,'/init']
save(name)

