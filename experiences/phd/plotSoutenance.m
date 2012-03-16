clear all

% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

set(0,'defaulttextinterpreter','tex')

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxeslinewidth',1);
set(0,'defaultlinelinewidth',2);

SizeMarker = 2;
SizeLabel = 30;

numfig = 1;
figure(numfig)

%%
%% Gibbs CirrusDot
%%

placemount = '/espace/orieux/results/'

load([placemount,'/init'])

expname = 'usmseCirrusDot'

alphaB = [10 555];
betaB = [10 570];

boundAlpha = [160 420];
boundBeta = [160 420];

boundBetal = [160 420];

name = [placemount,'variousRegCirrusDCT/gpac_',num2str(1)];
load(name, 'xchap')
dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

figplace = '/home/orieux/tex/planches/2009.11.16-soutenance-these@lss/figs/';
system(['mkdir -p ',figplace]);

%%

load([placemount,expname,'/EAP'])
sol = conv2(xEap(:,:,2),fgaussian,'same');

true = conv2(sky(:,:,2),fgaussian,'same');

load('/espace/orieux/results/fillCirrusDot/coadd','coaddCirrusDot')
load('/espace/orieux/results/fillCirrusDot/optim50000','fillmap')

coadded = conv2(fillmap(:,:,2),fgaussian,'same');

clf
imagesc(alpha(alphaB(1):alphaB(2)), beta(betaB(1):betaB(2)), ...
        sol(alphaB(1):alphaB(2),betaB(1):betaB(2)), dynamique)
axis image
axis xy; axis off
colormap(gray)

name = 'eapCirrusDot';
print(1,'-depsc',[figplace,'/',name])

clf
imagesc(alpha(alphaB(1):alphaB(2)), beta(betaB(1):betaB(2)), ...
        sky(alphaB(1):alphaB(2),betaB(1):betaB(2)), dynamique)
axis image
axis xy; axis off
colormap(gray)

name = 'cirrusDot';
print(1,'-depsc',[figplace,'/',name])

load('/espace/orieux/results/fillCirrusDot/coadd')

clf
imagesc(coaddCirrusDot(:,:,2)/Hrond360(1), dynamique)
axis image
axis xy; axis off
colormap(gray)

name = 'coadd';
print(1,'-depsc',[figplace,'/',name])

% l1CirrusDotEAP = sum2(abs(sky(:,:,2) - xEap(:,:,2)))/sum2(abs(sky(:,:,2))) 
% l2CirrusDotEAP = sum2((sky(:,:,2) - xEap(:,:,2)).^2)/sum2((sky(:,:,2).^2)) 

%%

name = [placemount,'fillCirrusDot/optim50000'];
load(name, 'fillmap')

coadded = conv2(fillmap(:,:,2),fgaussian,'same');

name = [placemount,'fillCirrusDot/coadd'];
coadded = coaddCirrusDot(:,:,2);

ligne = 340;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

%ynear = interp1 (xi, yi, x,'nearest');

clf
plot(beta, true(ligne,:), 'r')
hold on
plot(beta, sol(ligne,:))
plot(beta_mm, coadded(ligne_mm,:)/Hrond360(1),'k')
grid on
xlim([beta(betaB(1)) beta(betaB(2))])
ylim([4 13]*1e-5)

legend('vrai', 'SR', 'co-add')
legend('Location',[0.6 0.6 0.32 0.25])
%legend('boxoff')

grid on
xlabel('\beta')

name = 'sliceCirrusDotMoy';
print(1,'-depsc',[figplace,'/',name])

clf
plot(beta, true(ligne,:))
xlim([beta(betaB(1)) beta(betaB(2))])
ylim([4 13]*1e-5)
grid on
xlabel('\beta')

name = 'sliceCirrusDotMoyVrai';
print(1,'-depsc',[figplace,'/',name])

clf
plot(beta, sol(ligne,:))
xlim([beta(betaB(1)) beta(betaB(2))])
ylim([4 13]*1e-5)
grid on
xlabel('\beta')

name = 'sliceCirrusDotMoySR';
print(1,'-depsc',[figplace,'/',name])

clf
plot(beta_mm, coadded(ligne_mm,:)/Hrond360(1))
xlim([beta(betaB(1)) beta(betaB(2))])
ylim([4 13]*1e-5)
grid on
xlabel('\beta')

name = 'sliceCirrusDotMoyCoadd';
print(1,'-depsc',[figplace,'/',name])

%%
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
dspProposed = estimCircularPSD(ufft2(window.*xEap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspNoisefree = estimCircularPSD(ufft2(window.*noisefree(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);
dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)/Hrond360(1)),deltaF);

f1 = min(find(rf/2 >= 0.001));
f2 = min(find(rf/2 >= 0.15));

%%
clf
loglog(rf(f1:f2)/alpha_step*sqrt(2), dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step*sqrt(2), dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step*sqrt(2), dspFill(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step*sqrt(2), dspNoisefree(f1:f2),'--k')

xlim([0.005 0.11])
ylim([1e-17 1e-7])

grid on

xlabel('arcsec^{-1}')
ylabel('|S|^2')
legend('vrai', 'SR', 'co-add', 'conv')
legend('Location','SW')%[0.2 0.2 0.2 0.15])

loglog([0.044, 0.044], [1e-20 1e0])
loglog([0.044, 0.044]/2, [1e-20 1e0])
text(0.039,1e-8,'F_c')
text(0.019,1e-8,'F_e')

name = 'psdEapCirrusDot';
print(1,'-depsc',[figplace,'/',name])
%%

break

%%
%%

set(0,'defaultaxesfontsize',30);


iter = [1:length(gbChain(2,:))];
currentMean = cumsum(gbChain(2,:));

clf
loglog(gbChain(2,:))
hold on
loglog(gbChain(2,1),'x')
xlim([iter(1) iter(end)])
grid on

name = 'chaineGammaBCirrusDot';
print(1,'-depsc',[figplace,'/',name])

%%

iter = [400:400+length(squeeze(gbChain(2,400:end))) - 1];
gbChainC = squeeze(gbChain(2,400:end)); 
currentMeanC = cumsum(gbChainC)./[1:numel(gbChainC)];

clf
plot(iter,gbChainC)
hold on
plot(iter,currentMeanC,'r')
xlim([iter(1) iter(end)])
ylim([0.97 1.03]*1e6)

grid on

name = 'sousChaineGammaCCirrusDot';
print(1,'-depsc',[figplace,'/',name])

%%

clf
hist(gbChain(2,400:end),50)
xlim([0.97 1.03]*1e6)

name = 'histGammaBCirrusDot';
print(1,'-depsc',[figplace,'/',name])

meanGb = mean(gbChain(2,400:end))
stdGb = sqrt(mean(gbChain(2,400:end).^2) - meanGb^2)

%%
%%

iter = [1:length(squeeze(gxChain(1,2,:)))]';
currentMean = cumsum(squeeze(gxChain(1,2,:)));

clf
loglog(squeeze(gxChain(1,2,:)))
hold on
loglog(gxChain(1,2,1),'x')
xlim([iter(1) iter(end)])
grid on

name = 'chaineGammaXCirrusDot';
print(1,'-depsc',[figplace,'/',name])

%%

iter = [400:400+length(squeeze(gxChain(1,2,400:end))) - 1]';
gxChainC = squeeze(gxChain(1,2,400:end)); 
currentMeanC = cumsum(gxChainC)./[1:numel(gxChainC)]';

clf
plot(iter,gxChainC)
hold on
plot(iter,currentMeanC,'r')
ylim([2.35 2.65]*1e11)
xlim([iter(1) iter(end)])

grid on

name = 'sousChaineGammaXCirrusDot';
print(1,'-depsc',[figplace,'/',name])

%%

clf
hist(gxChain(1,2,400:end),50)
xlim([2.35 2.65]*1e11)

name = 'histGammaXCirrusDot';
print(1,'-depsc',[figplace,'/',name])

meanGx = mean(squeeze(gxChain(1,2,400:end)))
stdGx = sqrt(mean(gxChain(1,2,400:end).^2) - meanGx^2)

%%

close all
