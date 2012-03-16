clear all
% TODO
% - slice cirrus dot
% - psd cirrus dot

% format long
path(path,'../../libspire')
path(path,'../../simulator')
path(path,'../../utils')
path(path,'../../')
addpath(genpath('/home/orieux/codes/orieux_matlab_tb'))

figplace = '/home/orieux/tex/papers/aa/tex/figs2';
system(['mkdir -p ',figplace]);

set(0,'defaultaxesfontsize',25);
set(0,'defaulttextfontsize',25);
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);

SizeMarker = 2;
SizeLabel = 30;

numfig = 1;
figure(numfig)

%%
%% Cirrus

load('/home/space/phd_results/variousDirty/dirtymaps')
load('/home/space/phd_results/fillCirrus/optim50000','fillmap')

placemount = '/home/space/phd_results/';
expname = 'variousRegCirrusDCT';
Ntasks = 144;

alphaB = [5 555];
betaB = [7 570];
boundAlpha = [160 420];
boundBeta = [160 420];
boundBetal = [160 420];

name = [placemount,expname,'/gpac_',num2str(1)];
load(name, 'xchap')
dynamique = [min(min(xchap(:,:,2))) max(max(xchap(:,:,2)))];

load([placemount,expname,'/initialization'])

%%

clf
imagesc(sky(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

print('-depsc',[figplace,'/cirrus']);

%%
%% Comp madmap

clf
imagesc(madmap(:,:,2)/Hrond360(1),dynamique)
axis image; axis xy; axis off
colormap gray
colorbar

print('-depsc',[figplace,'/coaddCirrus']);

% 
%% EAP
placemount = '/home/space/phd_results/';
load([placemount,'/init']);
expname = 'usmseCirrus';
placemount = '/home/space/phd_results/';

load([placemount,expname,'/EAP'])
sol = conv2(xEap(:,:,2),fgaussian,'same');
true = conv2(sky(:,:,2),fgaussian,'same');
coadded = conv2(fillmap(:,:,2),fgaussian,'same');

% load('/home/space/phd_results/variousDirty/dirtymaps')
% %load('/home/space/phd_results/fillCirrus/coadd','coaddCirrus')
% load('/home/space/phd_results/fillCirrus/optim50000','fillmap')

clf
imagesc(alpha(alphaB(1):alphaB(2)), beta(betaB(1):betaB(2)), ...
        sol(alphaB(1):alphaB(2),betaB(1):betaB(2)), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

figplace = '/home/orieux/tex/papers/aa/tex/figs2';
print('-depsc',[figplace,'/eapCirrus'])

% 
%% Slice

ligne = 200;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

clf
plot(beta(betaB(1):betaB(2)), true(ligne,betaB(1):betaB(2)), 'r')
hold on
plot(beta(betaB(1):betaB(2)), sol(ligne,betaB(1):betaB(2)))
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
% plot(beta(betaB(1):betaB(2)), coadded(ligne,(betaB(1):betaB(2)))/Hrond360(1),'k')
grid on
xlim([beta(betaB(1)) beta(betaB(2))])
ylim(dynamique)

legend('True', 'Proposed', 'Co-add')
legend('Location',[0.6 0.6 0.32 0.25])

grid on

print('-depsc',[figplace,'/sliceCirrusMoy'])


%%
%% DSP

% Window Gaussienne
sigma = 30;

NalphaW = boundAlpha(2) - boundAlpha(1);
NbetaW = boundBeta(2) - boundBeta(1);

[ALPHA BETA] = ndgrid([-(NalphaW)/2:(NalphaW)/2],[-(NbetaW)/2:(NbetaW)/2]);
window = exp(-1/2 * (ALPHA.^2 + BETA.^2)/sigma^2);

% 
%% PSD

noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./Hrond360(1)));

deltaF = 0.004;
%% rf is reduced frequency
[dspTrue rf] = estimCircularPSD(ufft2(window.*sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspProposed = estimCircularPSD(ufft2(window.*xEap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspNoisefree = estimCircularPSD(ufft2(window.*noisefree(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);
dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)/Hrond360(1)),deltaF);

bestDspNormalNoisCirrus = dspProposed;

f1 = min(find(rf/2 >= 0.001));
f2 = min(find(rf/2 >= 0.15));

clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')

xlim([0.005 0.11])
ylim([1e-17 1e-7])

grid on

xlabel('f [1/arcsec]')
legend('True', 'Proposed', 'Co-add', 'Conv')
legend('Location',[0.2 0.2 0.32 0.25])

print('-depsc',[figplace,'/psdEapCirrus'])

%%
iter = [1:length(gbChain(2,:))];
currentMean = cumsum(gbChain(2,:));

clf
loglog(gbChain(2,:))
hold on
loglog(gbChain(2,1),'x')
xlim([iter(1) iter(end)])

grid on

print('-depsc',[figplace,'/chaineGammaB'])

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

print('-depsc',[figplace,'/sousChaineGammaB'])

%%

clf
hist(gbChain(2,400:end),50)
xlim([0.97 1.03]*1e6)

print('-depsc',[figplace,'/histGammaB'])

meanGb = mean(gbChain(2,400:end))
stdGb = sqrt(mean(gbChain(2,400:end).^2) - meanGb^2)

%%

iter = [1:length(squeeze(gxChain(1,2,:)))]';
currentMean = squeeze(cumsum(squeeze(gxChain(1,2,:))))./iter;

clf
loglog(squeeze(gxChain(1,2,:)))
hold on
loglog(gxChain(1,2,1),'x')
xlim([iter(1) iter(end)])
grid on

print('-depsc',[figplace,'/chaineGammaX'])

%%

iter = [400:400+length(squeeze(gxChain(1,2,400:end))) - 1]';
gxChainC = squeeze(gxChain(1,2,400:end)); 
currentMeanC = cumsum(gxChainC)./[1:numel(gxChainC)]';

clf
plot(iter,gxChainC)
hold on
plot(iter,currentMeanC,'r')
ylim([2.6 2.9]*1e11)
xlim([iter(1) iter(end)])

grid on

print('-depsc',[figplace,'/sousChaineGammaX'])

meanGx = mean(squeeze(gxChain(1,2,400:end)))
stdGx = sqrt(mean(gxChain(1,2,400:end).^2) - meanGx^2)

%%

clf
hist(gxChain(1,2,400:end),50)
xlim([2.6 2.9]*1e11)

print('-depsc',[figplace,'/histGammaX'])

%%
%% Variance
break

cumulant1 = zeros(size(sol));
cumulant2 = zeros(size(sol));

cumulant1r = zeros(size(dspProposed));
cumulant2r = zeros(size(cumulant1r));

hit = 0;
hist2D_DSP = zeros(100,100);
for itask = 400:932
    
    disp(num2str(itask));
    
    name = [placemount,expname,'/samples/sample_',num2str(itask)];
    try
        load(name, 'skySample')

        im = conv2(skySample(:,:,2),fgaussian,'same');
        imr = estimCircularPSD(ufft2(window.*im(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))), deltaF);
        
        cumulant1 = cumulant1 + im;
        cumulant1r = cumulant1r + imr;

        cumulant2 = cumulant2 + im.^2;
        cumulant2r = cumulant2r + imr.^2;
        
        [hist_dsp, bins] = hist3(log([imr(f1:f2), rf(f1:f2)/alpha_step]), [100 100]);
        hist2D_DSP = hist2D_DSP + hist_dsp;
        
        hit = hit+1;
    end
end

cumulant1 = cumulant1/hit;
cumulant2 = cumulant2/hit;

cumulant1r = cumulant1r/hit;
cumulant2r = cumulant2r/hit;

standartDev = sqrt(cumulant2 - cumulant1.^2);

standartDeviationR = sqrt(cumulant2r - cumulant1r.^2);

mean2(standartDev(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)))

clf
imagesc(alpha,beta,sqrt(cumulant2 - cumulant1.^2))
axis image
axis xy; axis off
colorbar
colormap(gray)

print('-depsc',[figplace,'/varianceCirrus'])

%%

clf
plot(beta(10:Nbeta), standartDev(ligne,10:Nbeta), 'k')
grid on
xlim([beta(10) beta(end)])

print('-depsc',[figplace,'/sliceStd'])


%%

clf
plot(beta(15:end-10), true(ligne,15:end-10))
hold on
plot(beta(15:end-10), sol(ligne,15:end-10) + 3*standartDev(ligne,15:end-10), 'k--')
plot(beta(15:end-10), sol(ligne,15:end-10) - 3*standartDev(ligne,15:end-10), 'k--')
grid on
xlim([beta(15) beta(end-10)])

print('-depsc',[figplace,'/sliceStandartDevGibbsCirrus'])

%%

clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, bestDspNormalNoisCirrus(f1:f2) + 3*standartDeviationR(f1:f2),'k--')
loglog(rf(f1:f2)/alpha_step, bestDspNormalNoisCirrus(f1:f2) - 3*standartDeviationR(f1:f2),'k--')
xlim([0.005 0.11])
ylim([1e-17 1e-7])
% xlim([1e-2 0.11])
% ylim([1e-14 1e-10])

grid on
xlabel('f [1/arcsec]')

print('-depsc',[figplace,'/psdCirrusStd']);

break

%% 
% Hist

clf
imagesc(bins{2}, bins{1}, log(hist2D_DSP))
hold on
plot(log(rf(f1:f2)/alpha_step), log(dspTrue(f1:f2)),'r')
axis xy
    
%%
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, bestDspNormalNoisCirrus(f1:f2) + 3*standartDeviationR(f1:f2),'k--')
loglog(rf(f1:f2)/alpha_step, bestDspNormalNoisCirrus(f1:f2) - 3*standartDeviationR(f1:f2),'k--')
xlim([0.005 0.11])
ylim([1e-17 1e-7])
% xlim([1e-2 0.11])
% ylim([1e-14 1e-10])

grid on
xlabel('f [1/arcsec]')

print('-depsc',[figplace,'/psdCirrusStd']);

break

% 
%% CirrusDot

load('/home/space/phd_results/fillCirrusDot/coadd','coaddCirrusDot')
madmap = coaddCirrusDot;
clear coaddCirrusDot;

boundAlpha = [160 420];
boundBeta = [160 420];

placemount = '/home/space/phd_results/';
expname = 'variousRegCirrusDotDCT';
Ntasks = 144;

load([placemount,expname,'/initialization'])

ligne = 340;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

%%

clf
imagesc(madmap(:,:,2)/Hrond360(1), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

print('-depsc',[figplace,'/coaddCirrusDot'])

%%

clf
imagesc(sky(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

print('-depsc',[figplace,'/cirrusDot'])

%%

placemount = '/home/space/phd_results/'
load([placemount,'/init'])
expname = 'usmseCirrusDot'

alphaB = [10 555];
betaB = [10 570];
boundAlpha = [160 420];
boundBeta = [160 420];
boundBetal = [160 420];

%%

load([placemount,expname,'/EAP'])
sol = conv2(xEap(:,:,2),fgaussian,'same');

true = conv2(sky(:,:,2),fgaussian,'same');

load('/home/space/phd_results/fillCirrusDot/coadd','coaddCirrusDot')
load('/home/space/phd_results/fillCirrusDot/optim50000','fillmap')

coadded = conv2(fillmap(:,:,2),fgaussian,'same');

clf
imagesc(alpha(alphaB(1):alphaB(2)), beta(betaB(1):betaB(2)), ...
        sol(alphaB(1):alphaB(2),betaB(1):betaB(2)), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)


print('-depsc',[figplace,'/eapCirrusDot'])

%%

name = [placemount,'fillCirrusDot/optim50000'];
load(name, 'fillmap')

coadded = conv2(fillmap(:,:,2),fgaussian,'same');

ligne = 340;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

clf
plot(beta(betaB(1):betaB(2)), true(ligne,betaB(1):betaB(2)), 'r')
hold on
plot(beta(betaB(1):betaB(2)), sol(ligne,betaB(1):betaB(2)))
plot(beta(betaB(1):betaB(2)), coadded(ligne,(betaB(1):betaB(2)))/Hrond360(1),'k')
grid on
ylim([4 13]*1e-5)

legend('True', 'Proposed', 'Co-add')
legend('Location',[0.6 0.6 0.32 0.25])

grid on

print('-depsc',[figplace,'/sliceCirrusDotMoy'])

%%

noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./Hrond360(1)));

deltaF = 0.004;
%% rf is reduced frequency
[dspTrue rf] = estimCircularPSD(ufft2(window.*sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspProposed = estimCircularPSD(ufft2(window.*xEap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspNoisefree = estimCircularPSD(ufft2(window.*noisefree(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);
dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)/Hrond360(1)),deltaF);

f1 = min(find(rf/2 >= 0.001));
f2 = min(find(rf/2 >= 0.15));

clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')

xlim([0.005 0.11])
ylim([1e-17 1e-7])

grid on

xlabel('f [1/arcsec]')
legend('True', 'Proposed', 'Co-add', 'Conv')
legend('Location',[0.2 0.2 0.32 0.25])

print('-depsc',[figplace,'/psdEapCirrusDot'])
