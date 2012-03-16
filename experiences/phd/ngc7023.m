clear all

% where result will be saved
placemount = '/espace/orieux/results/'
expname = 'ngc7023-2'
system(['mkdir -p ',placemount,expname]);

disp('Init')

% load data and pointing
load /espace/orieux/realData/data
load /espace/orieux/realData/pointing

pointing250 = pointing{1};
pointing360 = pointing{2};
pointing520 = pointing{3};

% set path, specially gpac
path(path,'/home/orieux/codes/spire/trunk')
path(path,'/home/orieux/codes/spire/trunk/libspire')
%path(path,'/home/orieux/codes/spire/trunk/simulator')
path(path,'/home/orieux/codes/spire/trunk/utils')
path(path,'/home/seismic/matlab/optimgpi')


%% Apply Alex formula to get a square map (abstraction of spherical
%coordinate ?)
for iscan = 1:10
    p = pointing250{iscan};
    p(1,:) = rad2arcsec(deg2rad(p(1,:)).*cos(7.5));
    p(2,:) = rad2arcsec(deg2rad(p(2,:)));
    pointing250(iscan) = {p};

    p = pointing360{iscan};
    p(1,:) = rad2arcsec(deg2rad(p(1,:)).*cos(7.5));
    p(2,:) = rad2arcsec(deg2rad(p(2,:)));
    pointing360(iscan) = {p};
    
    p = pointing520{iscan};
    p(1,:) = rad2arcsec(deg2rad(p(1,:)).*cos(7.5));
    p(2,:) = rad2arcsec(deg2rad(p(2,:)));
    pointing520(iscan) = {p};
end

% set parameters
paramsInstrument
paramsObservation
disp('The speed vector is not good !')
paramsSkyNGC7023

%% Precomputation of index and redondancy
[alpha, beta, pointing250, pointing360, pointing520, index250, index360, ...
 index520, coefs250, coefs360, coefs520] = computeIndex(alpha, beta, ...
                                                  pointing250, pointing360, ...
                                                  pointing520, alpha_step, ...
                                                  beta_step);

% Computation of coadd map. Need also new index and redondancy since the
% resolution is not the same.
disp('Dirty map')
[alpha_mm, beta_mm, pointing250_mm, pointing360_mm, pointing520_mm, ...
 index250_mm, index360_mm, index520_mm, coefs250_mm, coefs360_mm, coefs520_mm] ...
    = computeIndex(alpha, beta, pointing250, pointing360, pointing520, 12, 12);

% Compute for the three array
bands = [1 1 1];

% coaddition
coaddCirrus = dirtymap(data, bands, index250_mm, index360_mm, index520_mm, ...
                       coefs250_mm, coefs360_mm, coefs520_mm, ...
                       length(unique(alpha_mm)) , length(unique(beta_mm)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
%% Offsets
                   
%% Estimation of offset with only the pointing model.


offsets250 = zeros(1,Nbolo250);
offsets360 = zeros(1,Nbolo360);
offsets520 = zeros(1,Nbolo520);

offsets = {offsets250, offsets360, offsets520};

Nalpha_mm = length(unique(alpha_mm));
Nbeta_mm = length(unique(beta_mm));

disp('Start offset estimation')

for iter = 1:100

    coaddCirrus = dirtymapOffsets(data, bands, index250_mm, index360_mm, ...
                                  index520_mm, coefs250_mm, coefs360_mm, ...
                                  coefs520_mm, offsets, Nalpha_mm, Nbeta_mm);

    dataReproduction = directDirty(coaddCirrus, bands, index250_mm, ...
                                   index360_mm, index520_mm, Nalpha_mm, ...
                                   Nbeta_mm);

    offsets = estimOffsets(data, bands, dataReproduction);

end

disp('End offset estimation')

%% Correction of data

scans250 = data{1}; scans360 = data{2}; scans520 = data{3};

for iscan = 1:N_scan_total

    scans250(iscan) = {scans250{iscan} - repmat(offsets{1}, size(scans250{iscan},1), 1)};
    scans360(iscan) = {scans360{iscan} - repmat(offsets{2}, size(scans360{iscan},1), 1)};
    scans520(iscan) = {scans520{iscan} - repmat(offsets{3}, size(scans520{iscan},1), 1)};

end

data = {scans250, scans360, scans520};

%% END offsets.

%% New coadd
coaddCirrus2 = dirtymap(data, bands, index250_mm, index360_mm, index520_mm, ...
                        coefs250_mm, coefs360_mm, coefs520_mm, ...
                        length(unique(alpha_mm)) , length(unique(beta_mm)));

%% Coadd but with high resolution
trou = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                coefs360, coefs520, length(unique(alpha)), ...
                length(unique(beta)), '0');

%%

disp('Precalculs RI')
precalculsInvariantRI

%% Data simulation
std250 = (1e-3); gammaB250 = 1/std250^2;
std360 = std250; gammaB360 = 1/std360^2;
std520 = std360; gammaB520 = 1/std520^2;

std = [std250 std360 std520];
gammaB = [gammaB250, gammaB360, gammaB520];


%% GPAC Options
OPS(1)	 = 1;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 1e-3;	% Seuil d'arret sur x
OPS(3)	 = 1e-3;	% Seuil d'arret sur f
OPS(4)	 = 10;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 200;		% Nbre maximum d'itérations
OPS(15)	 = 1e-15;	% Pas minimum (MuEps)
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
hypers(2,:,1) = [1 1 1]*1e13;

%% Init
%%% The dirtymap filled with mean
init = dirtymap(data, [1 0 0], index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

%%% Gain adaptation
init(:,:,1) = init(:,:,1)./Hrond250(1);
init(:,:,2) = init(:,:,2)./Hrond360(1);
init(:,:,3) = init(:,:,3)./Hrond520(1);

%%% Prior mean
%%%% No prior mean here (only usefull for sampling)
objectMean = zeros(size(init));

%% Computation for 360 only
bands = [1 0 0];
bands = [1 1 1];

OPS(18)	 = 2e-12;		% Pas a l'origine

T0 = cputime;

disp('Opt')

% [xchap SortieOPS Histo] = gpac('calcQuadCrit', init, OPS,'calcQuadGrad', ...
%                                data, hypers, bands, Hrond250, Hrond360, ...
%                                Hrond520, index250, index360, index520, ...
%                                regOp, objectMean, Nalpha, Nbeta, Norder);

% DQ = cputime-T0;
% fprintf('\n\n\t Durée = %5.0f secondes i.e. %2.1f minutes \n\n',DQ,DQ/60)

% name = [placemount,expname,'/gpac0_360_',num2str(hypers(2,1,1))]
% save(name, 'xchap', 'OPS', 'SortieOPS', 'Histo', 'hypers', 'init')

%%

% Code for parallel computation settings

disp('Parallele')

Ntasks = 8*2;

variousReg = logspace(7,12,Ntasks);

save([placemount,expname,'/initialization'])

res = findResource('scheduler', 'type', 'jobmanager', 'name', 'jobmanager', ...
                   'LookupURL', 'cluster2.lss.supelec.fr');


%%PSW 

job = createJob(res);
% set path for job on cluster
set(job,'fileDependencies',{'.' '..' '../libspire' '../utils' ['/home/' ...
                    'seismic/matlab/optimgpi']})

disp('Init done')

disp('Submit tasks')

% build on task per regularisation parameters
for itask = 1:Ntasks
    % Diff for O1 only
    hypers(2,:,1) = variousReg(itask);

    bands = [1 0 0];
    
    theTask(itask) = createTask(job, 'gpac', 3, {'calcQuadCrit', init, ...
                        OPS,'calcQuadGrad', data, hypers, bands, Hrond250, ...
                        [], [], index250, [], [], regOp, objectMean, Nalpha, ...
                        Nbeta, Norder});

end

submit(job)

disp('submited')

waitForState(job,'finished')              

alltasks = get(job, 'Tasks');

% get results
for itask = 1:length(alltasks)
    
    result = get(alltasks(itask), 'OutputArguments');
    
    xchap = result{1};
    SortieOPS = result{2};
    Histo = result{3};
    hypers(2,:,1) = variousReg(itask);
    
    name = [placemount,expname,'/gpac_PSW',num2str(itask)];
    save(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

end

destroy(job)

%%PMW 

job = createJob(res);
set(job,'fileDependencies',{'.' '..' '../libspire' '../utils' ['/home/' ...
                    'seismic/matlab/optimgpi']})

disp('Init done')

disp('Submit tasks')

for itask = 1:Ntasks
    % Diff for O1 only
    hypers(2,:,1) = variousReg(itask);

    bands = [0 1 0];

    theTask(itask) = createTask(job, 'gpac', 3, {'calcQuadCrit', init, ...
                        OPS,'calcQuadGrad', data, hypers, bands, [], ...
                        Hrond360, [], [], index360, [], regOp, objectMean, ...
                        Nalpha, Nbeta, Norder});

end

submit(job)

disp('submited')

waitForState(job,'finished')              

alltasks = get(job, 'Tasks');

for itask = 1:length(alltasks)
    
    result = get(alltasks(itask), 'OutputArguments');
    
    xchap = result{1};
    SortieOPS = result{2};
    Histo = result{3};
    hypers(2,:,1) = variousReg(itask);
    
    name = [placemount,expname,'/gpac_PMW',num2str(itask)];
    save(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

end

destroy(job)

%%PLW 

job = createJob(res);
set(job,'fileDependencies',{'.' '..' '../libspire' '../utils' ['/home/' ...
                    'seismic/matlab/optimgpi']})

disp('Init done')

disp('Submit tasks')

for itask = 1:Ntasks
    % Diff for O1 only
    hypers(2,:,1) = variousReg(itask);

    bands = [0 0 1];

    theTask(itask) = createTask(job, 'gpac', 3, {'calcQuadCrit', init, ...
                        OPS,'calcQuadGrad', data, hypers, bands, [], [], ...
                        Hrond520, [], [], index520, regOp, objectMean, ...
                        Nalpha, Nbeta, Norder});

end

submit(job)

disp('submited')

waitForState(job,'finished')              

alltasks = get(job, 'Tasks');

for itask = 1:length(alltasks)
    
    result = get(alltasks(itask), 'OutputArguments');
    
    xchap = result{1};
    SortieOPS = result{2};
    Histo = result{3};
    hypers(2,:,1) = variousReg(itask);
    
    name = [placemount,expname,'/gpac_PLW',num2str(itask)];
    save(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

end

destroy(job)

clear res


