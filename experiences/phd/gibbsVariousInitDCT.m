%% variousRegCirruDotsDCT -  Different reg for cirrusDCT

placemount = '/espace/orieux/results/'
expname = 'gibbsVariousSkyDCT'
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
%OPS(9)(12)(13)(16)(17)(21)

% Gibbs options

criterion = 1e-4;
burnin = 400; %% If you absolutly don't knwo this value it is more than 100
maxIter = 1000; %% If you don't know it is more 200

%% Regularization

regOps = circDalpha + circDbeta;

%% Initial state of the chain


%% Initialisation (zero imply that the first conditionnal law is the hyper
%law. Non zero imply that the first conditionnal law is the image law)

hypersInit = zeros(2,3);
hypersInit(1,:) = 1.0e+06*1.0724;
hypersInit(2,:) = 2.5e11; %% Result from previous experience (it was
                          %convergend around this value)

%% computation for 360 only
bands = [0 1 0];

%%

%% DCT specific
res = findResource('scheduler', 'type', 'jobmanager', 'name', ...
                   'jobmanager', 'LookupURL', 'cluster2.lss.supelec.fr');

job = createJob(res);
set(job,'fileDependencies',{'.' '..' '../libspire' '../utils' '/home/seismic/matlab/optimgpi'})

%%

[data output data250 data360 data520 sky] = simulateData(['../data&true/' ...
                    'simul_ism_SPS_comp'], Nalpha, Nbeta, Norder, [1 1 1], ...
                                                  std, 10^(-4), Hrond250, ...
                                                  Hrond360, Hrond520, ...
                                                  index250, index360, ...
                                                  index520);

init = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

init = ones(size(init))*init(1,1,2);

%%% Adapte the gain. This is not a hack ! This is logical.
init(:,:,1) = init(:,:,1)./Hrond250(1);
init(:,:,2) = init(:,:,2)./Hrond360(1);
init(:,:,3) = init(:,:,3)./Hrond520(1);

theTask(1) = createTask(job, 'usmse', 3, {init, hypersInit, data, ...
                    bands, Hrond250, Hrond360, Hrond520, index250, index360, ...
                    index520, regOps,Nalpha, Nbeta, Norder, criterion, ...
                    burnin, maxIter, OPS});

%%

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
                                              Nalpha, Nbeta, Norder, [1 1 ...
                    1], std, 1, Hrond250, Hrond360, Hrond520, index250, ...
                                              index360, index520);

init = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

%%% Adapte the gain. This is not a hack ! This is logical.
init(:,:,1) = init(:,:,1)./Hrond250(1);
init(:,:,2) = init(:,:,2)./Hrond360(1);
init(:,:,3) = init(:,:,3)./Hrond520(1);

theTask(2) = createTask(job, 'usmse', 3, {init, hypersInit, data, ...
                    bands, Hrond250, Hrond360, Hrond520, index250, index360, ...
                    index520, regOps,Nalpha, Nbeta, Norder, criterion, ...
                    burnin, maxIter, OPS});

%%

[data output data250 data360 data520 sky] = simulateData('../data&true/simul_cond_ism_SPS_comp', ...
                                              Nalpha, Nbeta, Norder, [1 1 ...
                    1], std, 10^(-4), Hrond250, Hrond360, Hrond520, index250, ...
                                              index360, index520);

init = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

init = ones(size(init))*init(1,1,2);

%%% Adapte the gain. This is not a hack ! This is logical.
init(:,:,1) = init(:,:,1)./Hrond250(1);
init(:,:,2) = init(:,:,2)./Hrond360(1);
init(:,:,3) = init(:,:,3)./Hrond520(1);

theTask(3) = createTask(job, 'usmse', 3, {init, hypersInit, data, ...
                    bands, Hrond250, Hrond360, Hrond520, index250, index360, ...
                    index520, regOps,Nalpha, Nbeta, Norder, criterion, ...
                    burnin, maxIter, OPS});

%%

[data output data250 data360 data520 sky]] = simulateData('../data&true/galaxie', ...
                                              Nalpha, Nbeta, Norder, [1 1 ...
                    1], std, 10^(-4), Hrond250, Hrond360, Hrond520, index250, ...
                                              index360, index520);

init = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

init = zeros(size(init));

%%% Adapte the gain. This is not a hack ! This is logical.
init(:,:,1) = init(:,:,1)./Hrond250(1);
init(:,:,2) = init(:,:,2)./Hrond360(1);
init(:,:,3) = init(:,:,3)./Hrond520(1);

theTask(4) = createTask(job, 'usmse', 3, {init, hypersInit, data, ...
                    bands, Hrond250, Hrond360, Hrond520, index250, index360, ...
                    index520, regOps,Nalpha, Nbeta, Norder, criterion, ...
                    burnin, maxIter, OPS});

%%

submit(job)

disp('Submitted')

waitForState(job,'finished')              

alltasks = get(job, 'Tasks'); 

%%

result = get(alltasks(1), 'OutputArguments');

xEap = result{1};
gbChain = result{2};
gxChain = result{3};

name = [placemount,expname,'/cirrus'];
save(name, 'xEap', 'gbChain', 'gxChain')

%%

result = get(alltasks(2), 'OutputArguments');

xEap = result{1};
gbChain = result{2};
gxChain = result{3};

name = [placemount,expname,'/prior'];
save(name, 'xEap', 'gbChain', 'gxChain')

%%

result = get(alltasks(3), 'OutputArguments');

xEap = result{1};
gbChain = result{2};
gxChain = result{3};

name = [placemount,expname,'/cirrusDot'];
save(name, 'xEap', 'gbChain', 'gxChain')

%%

result = get(alltasks(4), 'OutputArguments');

xEap = result{1};
gbChain = result{2};
gxChain = result{3};

name = [placemount,expname,'/galaxie'];
save(name, 'xEap', 'gbChain', 'gxChain')

%%

destroy(job)

clear res



