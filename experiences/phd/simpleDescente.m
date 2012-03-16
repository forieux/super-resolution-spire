%% variousRegCirrusDCT - Compute several estimation with different reg

clear all

%% Base
% format long
path(path,'../../')
path(path,'../../libspire')
path(path,'../../data&true')
path(path,'../../utils')
path(path,'../../simulator')
path(path,'/Users/seismic/matlab/optimgpi')

%%

placemount = '/espace/orieux/results/'
expname = 'simpleCirrus'

%% 
paramsInstrument
paramsObservation; simulatePointing;
paramsSky

%% Precomputation
[alpha, beta, pointing250bis, index250, coefs250, Nalpha, Nbeta] = ...
        computeIndex(alpha, beta, pointing250, alpha_step, beta_step, ...
                      Nbolo250, N_scan_total);


[Hrond250 Hdirect250] = computeRI(10, 1, band_250, ...
                                  central_wavelength_250, ...
                                  sigma_coef, sigma_alpha, ...
                                  sigma_beta, alpha_step, beta_step, ...
                                  Nalpha, Nbeta, Nspeed, Norder, ...
                                  unique_speed);

[diffAlpha diffBeta circDalpha circDbeta] = ...
      computeReg(sigma_alpha, sigma_beta, sky_alpha_period, ...
                 sky_beta_period, 5, Nalpha, Nbeta);
               
%% Data simulation
std250 = (1e-3); gammaB250 = 1/std250^2;
std360 = std250; gammaB360 = 1/std360^2;
std520 = std360; gammaB520 = 1/std520^2;

std = [std250 std360 std520];
gammaB = [gammaB250, gammaB360, gammaB520];

[data250 output250 sky] = simulateData('../../data&true/simul_ism_SPS_comp', Nalpha, Nbeta, Norder, ...
                                                  std250, 10^(-4), Hrond250, ...
                                                  index250, Nbolo250, Nspeed, N_scan_total, unique_speed, the_speeds);

%% GPAC Options
OPS(1)	 = 20;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 1e-6;	% Seuil d'arret sur x
OPS(3)	 = 1e-6;	% Seuil d'arret sur f
OPS(4)	 = 10000;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 50;		% Nbre maximum d'itérations
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
hypers(2,:,1) = 1.4*1e12;


%% Init
%%% The dirtymap filled with mean
init = zeros(Nalpha, Nbeta, 3);

% init = dirtymap(data, index, index360, index520, coefs250, ...
%                  coefs360, coefs520, Nalpha, Nbeta, 'M');

%%% Or load previous computed filled map (see simpleFill)
% load /espace/orieux/results/simpleFill/optim500000 fillmap;
% init = fillmap;

%%% Gain adaptation
% init(:,:,1) = init(:,:,1)./Hrond250(1);
% init(:,:,2) = init(:,:,2)./Hrond360(1);
% init(:,:,3) = init(:,:,3)./Hrond520(1);

%%% Prior mean
%%%% No prior mean here (only usefull for sampling)
init = zeros(Nalpha, Nbeta);
objectMean = zeros(size(init));

%% Computation for 360 only
bands = [0 1 0];

OPS(18)	 = 1e-12;		% Pas a l'origine

tic
[xchapGPAC SortieOPS Histo] = gpac('calcQuadCrit', init, OPS,'calcQuadGrad', ...
                               data, hypers, bands, Hrond250, Hrond360, ...
                               Hrond520, index250, index360, index520, ...
                               regOp, objectMean, Nalpha, Nbeta, Norder);
TGPAC = toc;
critGPAC = calcQuadCrit(xchapGPAC, data, hypers, bands, Hrond250, Hrond360, ...
                             Hrond520,  index250, index360, index520, ...
                             regOp, objectMean, Nalpha, Nbeta, Norder);

dataProj = transposeInvariant(data, bands, Hrond250, Hrond360, Hrond520, ...
                              index250, index360, index520, Nalpha, Nbeta, ...
                              Norder);
dataProj = myfft2(dataProj);
init = myfft2(init);

dataProj(:,:,1:3:end) = 2*hypers(1,1,1)*dataProj(:,:,1:3:end);
dataProj(:,:,2:3:end) = 2*hypers(1,2,1)*dataProj(:,:,2:3:end);
dataProj(:,:,3:3:end) = 2*hypers(1,3,1)*dataProj(:,:,3:3:end);

cgoptions.thresold = 1e-1000; cgoptions.maxIter = OPS(14); %cgoptions.numfig = 1000;

scoefs250 = zeros(Nalpha, Nbeta, Nspeed);
scoefs360 = zeros(Nalpha, Nbeta, Nspeed);
scoefs520 = zeros(Nalpha, Nbeta, Nspeed);

for iscan = 1:N_scan_total
  
  %% identify the PSF with the speed
  speed_index = find(unique_speed(1,:) == the_speeds(1,1,iscan));
  
  scoefs250(:,:,speed_index) = scoefs250(:,:,speed_index) + ...
      coefs250{iscan};
  
  scoefs360(:,:,speed_index) = scoefs360(:,:,speed_index) + ...
      coefs360{iscan};
  
  scoefs520(:,:,speed_index) = scoefs520(:,:,speed_index) + ...
      coefs520{iscan};
  
end

tic
[xchapCG SortieOPS Histo] = conjGrad('appHessian', init, dataProj , ...
                                   cgoptions, 'precond', hypers, ...
                                   bands, Hrond250, Hrond360, ...
                                   Hrond520, scoefs250, scoefs360, ...
                                   scoefs520, regOp, objectMean, ...
                                   Nalpha, Nbeta, Norder);

TCG = toc;
critCG = calcQuadCrit(real(myifft2(xchapCG)), data, hypers, bands, Hrond250, Hrond360, ...
                             Hrond520,  index250, index360, index520, ...
                             regOp, objectMean, Nalpha, Nbeta, Norder);

disp(['Temps GPAC : ',num2str(TGPAC)])
disp(['Temps CG : ',num2str(TCG)])
disp(['Crit GPAC : ',num2str(critGPAC)])
disp(['Crit CG : ',num2str(critCG)])
