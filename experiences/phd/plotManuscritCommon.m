clear all

% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

figplace = '/home/orieux/tex/manuscrit/imagesRes/common/';
system(['mkdir -p ',figplace,'/eps']);

set(0,'defaulttextinterpreter','none')

set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxeslinewidth',1);
set(0,'defaultlinelinewidth',0.5);

SizeMarker = 2;
SizeLabel = 30;

numfig = 1;
figure(numfig)

%%
%% Cirrus
%%

load('/espace/orieux/results/variousDirty/dirtymaps')
load('/espace/orieux/results/fillCirrus/optim50000','fillmap')

placemount = '/espace/orieux/results/'
expname = 'commonSig';
Ntasks = 144;

boundAlpha = [160 420];
boundBeta = [160 420];

boundBetal = [160 420];


name = [placemount,'variousRegCirrusDCT','/gpac_',num2str(1)];
load(name, 'xchap')
dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

load([placemount,'/commonSigResults'])

load([placemount,'init'], 'Hrond360')

[theImage comSig] = unpackSkyCommon(xchap, Nalpha, Nbeta, Norder);

set(0,'defaulttextinterpreter','none')

%%

clf
imagesc(theImage(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

name = 'cirrusCommon';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

ligne = 200;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), theImage(ligne,10:Nbeta,2))
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)
legend('vrai', 'mcr')
legend('Location',[0.6 0.6 0.32 0.25])

name = 'sliceCirrusCommon';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
imagesc(nocorrel(:,:,2)./max2(nocorrel(:,:,2))*max2(theImage(:,:,2)), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

name = 'cirrusCommonPaf';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

ligne = 200;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), nocorrel(ligne,10:Nbeta,2)./max2(nocorrel(:,:,2))*max2(theImage(:,:,2)))
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)
legend('vrai', 'mcr')
legend('Location',[0.6 0.6 0.32 0.25])

name = 'sliceCirrusCommonPaf';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
plot([1:1500].*temporal_sampling_periode,commonSignal,'r')
hold on
plot([1:1500].*temporal_sampling_periode,comSig(:,2))
xlim([1 1500].*temporal_sampling_periode)
ylim([0.2 0.65])
grid on
xlabel('$t~[\textrm{s}]$')
ylabel('$[\textrm{V}]$')

name = 'commonSignal';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

break

%%
%% DSP

%%
% Window

%% Gaussienne
sigma = 30;

NalphaW = boundAlpha(2) - boundAlpha(1);
NbetaW = boundBeta(2) - boundBeta(1);

[ALPHA BETA] = ndgrid([-(NalphaW)/2:(NalphaW)/2],[-(NbetaW)/2:(NbetaW)/2]);
window = exp(-1/2 * (ALPHA.^2 + BETA.^2)/sigma^2);

%% cosinus

% Nbord = 30;

% bordw = cos(linspace(-pi/2,pi/2,Nbord));

% windowA = [bordw(1:Nbord/2) ones(NalphaW - Nbord + 1,1)' bordw(Nbord/2+1:end)];
% windowB = [bordw(1:Nbord/2) ones(NbetaW - Nbord + 1,1)' bordw(Nbord/2+1:end)];

%window = kron(windowA',windowB);
%window = ones(size(window));

% imagesc(window.*sky(boundAlpha(1):boundAlpha(2),boundBeta(1): boundBeta(2),2))

%

noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./Hrond360(1)));

deltaF = 0.004;
%% rf is reduced frequency
[dspTrue rf] = estimCircularPSD(ufft2(window.*sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspProposed = estimCircularPSD(ufft2(window.*image(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspNoisefree = estimCircularPSD(ufft2(window.*noisefree(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);
dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)/Hrond360(1)),deltaF);

f1 = min(find(rf/2 >= 0.001));
f2 = min(find(rf/2 >= 0.15));

clf
% RF = F/FE = F*Te

loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')
xlim([0 rf(f2)/alpha_step])
ylim([1e-25 1e-5])
grid on

xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
ylabel('$|S|^2$')
legend('vrai', 'mcr', 'co-add', 'conv')
legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')

title(['$\gx = 10^7$'])

name = 'psdLowRegCirrus';
laprint(1,[figplace,'/',name], 'width',8, 'asonscreen', 'on')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

