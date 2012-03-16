% Time-stamp: < 15/04(avr)/2010 15:00 by orieux (simpleWienerNS.m) >

% clear all

%% 

placemount = '/espace/orieux/results/'
expname = 'simpleWienerNS'
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

%% Precomputation
[alpha, beta, start_end_position, pointing250, pointing360, pointing520, ...
 index250, index360, index520, coefs250, coefs360, coefs520] = ...
    computeIndex(alpha, beta, start_end_position, pointing250, pointing360, ...
                 pointing520, alpha_step, beta_step);

precalculsInvariantRI

%% Data simulation
std250 = (1e-3); gammaB250 = 1/std250^2;
std360 = std250; gammaB360 = 1/std360^2;
std520 = std360; gammaB520 = 1/std520^2;

std = [std250 std360 std520];
gammaB = [gammaB250, gammaB360, gammaB520];

[data output data250 data360 data520 sky] = simulateData(['../data&true/' ...
                    'simul_ism_SPS_comp'], Nalpha, Nbeta, Norder, [1 1 1], ...
                                                  std, 10^(-4), Hrond250, ...
                                                  Hrond360, Hrond520, ...
                                                  index250, index360, ...
                                                  index520);

% hyper -- the hyperparamter value tab of 3*Norder + 1 lines and 3 columns,
% one column for each band. The line are ordered like this : the first line
% is noise precision. The lines 2 to 2+(Norder-1) is for difference in
% alpha. From 3+(Norder-1) to 3+2*(Norder-1) is for difference in beta. From
% 4+2*Norder to 4+3*(Norder-1) lines is for the mean. From 5+3*(Norder-1) to
% end is for the norm.

%%% Only on operator in line an column in same time
regOp = circDalpha + circDbeta;

%% Bound
alphaBound = zeros(2,3);
betaBound = zeros(2,3);

%% Data retroprojection
cielRondSmooth = zeros(Nalpha,Nbeta,3);

load /espace/orieux/results/simpleFillCirrus/optim500000 fillmap;

for iarray = 1:3
    cielRondSmooth(:,:,iarray) = ufft2(fillmap(:,:,iarray));
end

%% Estimation
bands = [0 1 0];
[xchap gbChain gChain] = wiener0NS(cielRondSmooth, bands, Hrond250, Hrond520, ...
                                   Hrond360, regOp, Nalpha, Nbeta, alphaBound, ...
                                   betaBound, 5e-4, 300, 1000, 'fig', 1000);

expname = 'simpleWienerNS';
name = [placemount,expname,'/EAP']
system(['mkdir -p ',placemount,expname]);

save(name, 'xchap', 'gbChain', 'gChain')

mgb = mean(gbChain,2)
mgx = mean(gChain(1,2,:))

lambda = mgx./mgb

