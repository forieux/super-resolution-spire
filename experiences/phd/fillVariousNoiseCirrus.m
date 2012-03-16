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
expname = 'FillVariousNoiseCirrus'
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

precalculsInvariantRI

%% GPAC Options
OPS(1)	 = 1;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 1e-6;	% Seuil d'arret sur x
OPS(3)	 = 1e-6;	% Seuil d'arret sur f
OPS(4)	 = 10;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 300;		% Nbre maximum d'itérations
OPS(15)	 = 1e-20;	% Pas minimum (MuEps)
OPS(18)	 = 1e-12;		% Pas a l'origine
%OPS(19) = 10;		% Sauve toutes le OPS(19) itérations...
%OPS(20) = 1;		% ... dans un fichier portant un numéro
%OPS(22) = 1;		% Positivité

%% Hyperparameters value
% hyper -- the hyperparamter value tab of 3*Norder + 1 lines and 3 columns,
% one column for each band. The line are ordered like this : the first line
% is noise precision. The lines 2 to 2+(Norder-1) is for difference in
% alpha. From 3+(Norder-1) to 3+2*(Norder-1) is for difference in beta. From
% 4+2*Norder to 4+3*(Norder-1) lines is for the mean. From 5+3*(Norder-1) to
% end is for the norm.

%%% Only on operator in line an column in same time
regOp = circDalpha + circDbeta;

%% Calcul 300 échantillons en //

variousNoise = [1e-2 1e-4];

save([placemount,expname,'/initialization'])

% %%

%% Data simulation
std250 = variousNoise(1); gammaB250 = 1/std250^2;
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

%% Regularization
hypers = zeros(2,3);
% Noise
hypers(1,:) = [gammaB250 gammaB360 gammaB520];
hypers(1,:) = [gammaB250 gammaB360 gammaB520];
% Diff (for [-1 1]) (5e7 is better for reg on gaussian)
hypers(2,1) = 3e4;
hypers(2,2) = 3e4;
hypers(2,3) = 3e4;

regOp = diffAlpha + diffBeta;

%% Init
init = dirtymap(data, [1 1 1], index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

bands = [0 1 0];

[fillmap SortieOPS Histo] = gpac('calcDirtyCrit',init,OPS,'calcDirtyGrad', ...
                                 data, hypers, bands, index250, index360, ...
                                 index520, regOp, Nalpha, Nbeta, 'fig', 1000);

name = [placemount,expname,'/highNoise',num2str(hypers(2,1))]
save(name, 'fillmap', 'OPS', 'SortieOPS', 'Histo', 'hypers', 'init')

%

%% Data simulation
std250 = variousNoise(2); gammaB250 = 1/std250^2;
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

%% Regularization
hypers = zeros(2,3);
% Noise
hypers(1,:) = [gammaB250 gammaB360 gammaB520];
hypers(1,:) = [gammaB250 gammaB360 gammaB520];
% Diff (for [-1 1]) (5e7 is better for reg on gaussian)
hypers(2,1) = 6e2;
hypers(2,2) = 6e2;
hypers(2,3) = 6e2;

regOp = diffAlpha + diffBeta;

%% Init
init = dirtymap(data, [1 1 1], index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

bands = [0 1 0];

[fillmap SortieOPS Histo] = gpac('calcDirtyCrit',init,OPS,'calcDirtyGrad', ...
                                 data, hypers, bands, index250, index360, ...
                                 index520, regOp, Nalpha, Nbeta, 'fig', 1000);

name = [placemount,expname,'/lowNoise',num2str(hypers(2,1))]
save(name, 'fillmap', 'OPS', 'SortieOPS', 'Histo', 'hypers', 'init')
