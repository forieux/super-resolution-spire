%% VARIANCECIRRUSDCT - Compute posterior variance of cirrus 

clear all

%% 

placemount = '/espace/orieux/results/'
expname = 'optimizeC'
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
%

%% Correlated signal simulation

%%% Same seed (see >> doc randn)
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));

commonSignal = randn(250*N_scan_total,1);

%%% DSP
frequency = linspace(0,1,length(commonSignal)); % reduced frequency
psdCommon = psdCommonSig(frequency, frequencyCutE, inclination);

%%% Filter in fourier
commonSignal = stdC*real(myifft(myfft(commonSignal').*sqrt(psdCommon)));
%%%% It is positive
commonSignal = commonSignal + 2*abs(min(commonSignal));

%%% ComputeIndex of commonSignal that match the used data used to estimate
computeDataIndex

data250 = output{1};
data360 = output{2};
data520 = output{3};

output250 = output{1};
output360 = output{2};
output520 = output{3};

for iscan = 1:N_scan_total
    common250_for_these_scan = repmat(commonSignal(indexCommon250(iscan,:)), ...
                                      [], Nbolo250);
    
    data250(iscan) = {data250{iscan} + common250_for_these_scan' + ...
                      std250*randn(size(data250{iscan}))};

    common360_for_these_scan = repmat(commonSignal(indexCommon360(iscan,:)), ...
                                      [], Nbolo360);

    data360(iscan) = {data360{iscan} + common360_for_these_scan' + ...
                      std360*randn(size(data360{iscan}))};

    common520_for_these_scan = repmat(commonSignal(indexCommon520(iscan,:)), ...
                                      [], Nbolo520);
        
    data520(iscan) = {data520{iscan} + common520_for_these_scan' + ...
                      std520*randn(size(data520{iscan}))};
end

for iscan = 1:N_scan_total
    data250NoCorrel(iscan) = {output250{iscan} + std250* ...
                        randn(size(data250{iscan}))};

    data360NoCorrel(iscan) = {output360{iscan} + std360* ...
                        randn(size(data360{iscan}))};

    data520NoCorrel(iscan) = {output520{iscan} + std520* ...
                        randn(size(data520{iscan}))};
end

data = {data250 data360 data520};
dataNoCorrel = {data250NoCorrel data360NoCorrel data520NoCorrel};


%% Blind Bolometer
dataBlindA250 = commonSignal + std250*randn(size(commonSignal));
dataBlindB250 = commonSignal + std250*randn(size(commonSignal));
             
dataBlindA360 = commonSignal + std360*randn(size(commonSignal));
dataBlindB360 = commonSignal + std360*randn(size(commonSignal));
             
dataBlindA520 = commonSignal + std520*randn(size(commonSignal));
dataBlindB520 = commonSignal + std520*randn(size(commonSignal));

dataBlind = {dataBlindA250 dataBlindB250 dataBlindA360 dataBlindB360 ...
             dataBlindA520 dataBlindB520};

%
%% GPAC Options
OPS(1)	 = 1;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 5e-10;	% Seuil d'arret sur x
OPS(3)	 = 5e-10;	% Seuil d'arret sur f
OPS(4)	 = 10;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 200;		% Nbre maximum d'itérations
OPS(15)	 = 1e-25;	% Pas minimum (MuEps)
OPS(18)	 = 1e-12;		% Pas a l'origine
%OPS(19) = 10;		% Sauve toutes le OPS(19) itérations...
%OPS(20) = 1;		% ... dans un fichier portant un numéro
%OPS(22) = 1;		% Positivité

initC = zeros(length(commonSignal),3);

initC(:,1) = mean([dataBlindA250;dataBlindB250]);
initC(:,2) = mean([dataBlindA360;dataBlindB360]);
initC(:,3) = mean([dataBlindA520;dataBlindB520]);

meanC = zeros(size(initC));
meanC(:,2) = randn(size(meanC(:,2)));
meanC(:,2) = stdC*real(myifft(myfft(meanC(:,2)').*sqrt(psdCommon)));

hypersC = zeros(2,3);
hypersC(1,:) = [gammaB250 gammaB360 gammaB520];
% hypersC(2,:) = 0;
% hypersC(2,:) = 1e-10;
hypersC(2,:) = 1/stdC^2;

OPS(18)	 = 1e-4;

bands = [0 1 0];

[xchap SortieOPS Histo] = gpac('calcCommonCrit', initC, OPS,'calcCommonGrad', ...
                               output, data, dataBlind, hypersC, bands, ...
                               indexCommon, psdCommon', meanC, 'fig', 1000);
figure(1)
clf
plot(xchap(:,2))
hold on
plot(commonSignal,'r')

% name = [placemount,expname,'/gpac0_360_',num2str(hypers(2,1,1))]
% save(name, 'xchap', 'OPS', 'SortieOPS', 'Histo', 'hypers', 'init')

