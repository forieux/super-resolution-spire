clear all

%% 

placemount = '/espace/orieux/results/'
expname = 'fillPrior'
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

%%% Only on operator in line an column in same time
regOp = circDalpha + circDbeta;

%%% Prior sample
gammaX = 4e11;

%%% Same seed (see >> doc randn)
RandStream.setDefaultStream(RandStream('mt19937ar','seed',0));

sky = uifft2(ufft2(randn(Nalpha,Nbeta))./(sqrt(gammaX*regOp)));
%%%% positivity
sky = sky + abs(min2(sky));

%%%% Save for simulateData
save('../data&true/prior','sky')
save('../data&true/gammaXprior','gammaX')

[data output data250 data360 data520 sky] = simulateData('../data&true/prior', ...
                                                  Nalpha, Nbeta, Norder, ...
                                                  [1 1 1], std, 1, Hrond250, ...
                                                  Hrond360, Hrond520, ...
                                                  index250, index360, ...
                                                  index520);
%% GPAC Options
OPS(1)	 = 1;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 2e-5;	% Seuil d'arret sur x
OPS(3)	 = 2e-5;	% Seuil d'arret sur f
OPS(4)	 = 10;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 300;		% Nbre maximum d'itérations
OPS(15)	 = 1e-15;	% Pas minimum (MuEps)
OPS(18)	 = 1e-10;	% Pas a l'origine
%OPS(19) = 10;		% Sauve toutes le OPS(19) itérations...
%OPS(20) = 1;		% ... dans un fichier portant un numéro

%% Regularization
hypers = zeros(2,3);
% Noise
hypers(1,:) = [gammaB250 gammaB360 gammaB520];
hypers(1,:) = [gammaB250 gammaB360 gammaB520];
% Diff (for [-1 1]) (5e7 is better for reg on gaussian)
hypers(2,1) = 5e6;
hypers(2,2) = 5e6;
hypers(2,3) = 5e6;

regOp = diffAlpha + diffBeta;

%% Init
init = dirtymap(data, [1 1 1], index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

bands = [1 1 1];

T0 = cputime;
[fillmap SortieOPS Histo] = gpac('calcDirtyCrit',init,OPS,'calcDirtyGrad', ...
                                 data, hypers, bands, index250, index360, ...
                                 index520, regOp, Nalpha, Nbeta, 'fig', 1000);
DQ = cputime-T0;
fprintf('\n\n\t Durée = %5.0f secondes i.e. %2.1f minutes \n\n',DQ,DQ/60)

name = [placemount,expname,'/optim',num2str(hypers(2,1))]
save(name, 'fillmap', 'OPS', 'SortieOPS', 'Histo', 'hypers', 'init')

%%

[alpha_mm, beta_mm, start_end_position_mm, pointing250_mm, pointing360_mm, ...
 pointing520_mm, index250_mm, index360_mm, index520_mm, coefs250_mm, ...
 coefs360_mm, coefs520_mm] = computeIndex(alpha, beta, start_end_position, ...
                                          pointing250, pointing360, ...
                                          pointing520, 6, 6);

coaddPrior = dirtymap(data, bands, index250_mm, index360_mm, index520_mm, ...
                      coefs250_mm, coefs360_mm, coefs520_mm, ...
                      length(unique(alpha_mm)) , length(unique(beta_mm)));

name = [placemount,expname,'/coadd']
save(name, 'coaddPrior', 'alpha_mm', 'beta_mm')
