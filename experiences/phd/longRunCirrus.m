%% VARIANCECIRRUSDCT - Compute posterior variance of cirrus 

clear all

%% 

placemount = '/espace/orieux/results/'
expname = 'longRunCirrus'
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
%% GPAC Options
OPS(1)	 = 1;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 0;	% Seuil d'arret sur x
OPS(3)	 = 0;	% Seuil d'arret sur f
OPS(4)	 = 10;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 200000;		% Nbre maximum d'itérations
OPS(15)	 = 0;	% Pas minimum (MuEps)
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

% Noise
hypers = zeros(2,3,1);
hypers(1,:,1) = [gammaB250 gammaB360 gammaB520];
% Diff for O1 only
hypers(2,:,1) = 3.9e12;

%% Init
%%% The dirtymap filled with mean
init = dirtymap(data, [1 1 1], index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

%%% Or load previous computed filled map (see simpleFill)
% load /espace/orieux/results/simpleFill/optim500000 fillmap;
% init = fillmap;

%%% Gain adaptation
init(:,:,1) = init(:,:,1)./Hrond250(1);
init(:,:,2) = init(:,:,2)./Hrond360(1);
init(:,:,3) = init(:,:,3)./Hrond520(1);

%%% Prior mean
%%%% No prior mean here (only usefull for sampling)
objectMean = zeros(size(init));

%% Computation for 360 only
bands = [0 1 0];

save([placemount,expname,'/initialization'])

[xchap SortieOPS Histo] = gpac('calcQuadCrit', init, OPS,'calcQuadGrad', ...
                               data, hypers, bands, Hrond250, Hrond360, ...
                               Hrond520, index250, index360, index520, ...
                               regOp, objectMean, Nalpha, Nbeta, Norder);

name = [placemount,expname,'/gpac']
save(name, 'xchap', 'OPS', 'SortieOPS', 'Histo', 'hypers', 'init')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% phi = conv2(xchap(:,:,2),fgaussian,'same');
% figure
% plot(phi(200,:))
% hold on
% plot(sky(200,:,2),'r')
% plot(mapmean(200,:,2)/Hrond360(1),'k')
