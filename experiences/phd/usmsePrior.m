%% variousRegCirruDotsDCT -  Different reg for cirrusDCT

placemount = '/espace/orieux/results/'
expname = 'usmsePrior'
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
OPS(1)	 = 5;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 1e-5;	% Seuil d'arret sur x
OPS(3)	 = 1e-5;	% Seuil d'arret sur f
OPS(4)	 = 7;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 200;		% Nbre maximum d'itérations
OPS(15)	 = 1e-20;	% Pas minimum (MuEps)
OPS(18)	 = 1e-10;		% Pas a l'origine
%OPS(19) = 10;		% Sauve toutes le OPS(19) itérations...
%OPS(20) = 1;		% ... dans un fichier portant un numéro
%OPS(22) = 1;		% Positivité

% Paramètres en sortie...
%OPS(8)  = Valeur finale du critère
%OPS(10) = Nbre d'évaluation du critère
%OPS(11) = Nbre d'évaluation du gradient
% Paramètres inutilisé
%OPS(9)(12)(13)(16s)(17)(21)

% Gibbs options

criterion = 1e-4;
burnin = 400; %% If you absolutly don't knwo this value it is more than 100
maxIter = 1000; %% If you don't know it is more 200

%% Regularization

regOps = circDalpha + circDbeta;

%% Initial state of the chain

%%% The dirtymap filled with mean
% init = dirtymap(data, bands, index250, index360, index520, coefs250, ...
%                      coefs360, coefs520, Nalpha, Nbeta, 'm');

%%% Or load previous computed filled map (see simpleFill)
load /espace/orieux/results/fillPrior/optim5000000 fillmap;
init = fillmap;
 
%%% Adapte the gain. This is not a hack ! This is logical.
init(:,:,1) = init(:,:,1)./Hrond250(1);
init(:,:,2) = init(:,:,2)./Hrond360(1);
init(:,:,3) = init(:,:,3)./Hrond520(1);

%% Initialisation (zero imply that the first conditionnal law is the hyper
%law. Non zero imply that the first conditionnal law is the image law)

hypersInit = zeros(2,3);
% hypersInit(1,:) = 1.0e+06*[1.0724, 0.6575, 0.2622]';
% hypersInit(2,:) = 2.78e11; %% Result from previous experience (it was
%                           %convergend around this value)
%hypersInit(2,:) = 1.0e+06*[0.2107, 1.1735, 0.0375]';

%% computation for 360 only
bands = [0 1 0];

[xEap gbChain gxChain] = usmse(init, hypersInit, data, bands, Hrond250, ...
                               Hrond360, Hrond520, index250, index360, ...
                               index520, regOps, Nalpha, Nbeta, Norder, ...
                               criterion, burnin, maxIter, OPS, 'pla', ...
                               [placemount, expname], 'fig', 1000);

name = [placemount,expname,'/EAP'];
try
    save(name, 'sky', 'init', 'hypersInit', 'xEap', 'gbChain', 'gxChain')
end

% %%% Certainly rubbish
% % %% Compute reg only with sky inside this bound
% b = bound/2;
% [I J] = ind2sub([Nalpha Nbeta], find(isnan(mapnan(:,:,1)) == 0));
% alphaBound(:,1) = [min(I)-b max(I)+b];
% betaBound(:,1) = [min(J)-b max(J)+b];

% [I J] = ind2sub([Nalpha Nbeta], find(isnan(mapnan(:,:,2)) == 0));
% alphaBound(:,2) = [min(I)-b max(I)+b];
% betaBound(:,2) = [min(J)-b max(J)+b];

% [I J] = ind2sub([Nalpha Nbeta], find(isnan(mapnan(:,:,3)) == 0));
% alphaBound(:,3) = [min(I)-b max(I)+b];
% betaBound(:,3) = [min(J)-b max(J)+b];

