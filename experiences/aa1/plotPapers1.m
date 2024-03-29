clear all

% format long
path(path,'../../libspire')
path(path,'../../simulator')
path(path,'../../utils')
path(path,'../../')
%path(path,'/home/seismic/matlab/optimgpi')

figplace = '/home/orieux/tex/papers/aa/tex/figs1';
system(['mkdir -p ',figplace]);

set(0,'defaultaxesfontsize',25);
set(0,'defaulttextfontsize',25);
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);

SizeMarker = 2;
SizeLabel = 30;

numfig = 1;
figure(numfig)

% imagesc(alpha, beta, sky, ); setim

% print('-depsc',[figplace,'/galaxieSat']);


%%
%% Cirrus

load('/media/espace/orieux/results/variousDirty/dirtymaps')
load('/media/espace/orieux/results/fillCirrus/optim50000','fillmap')

placemount = '/media/espace/orieux/results/';
expname = 'variousRegCirrusDCT';
Ntasks = 144;

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
%% L2 et L1
placemount = '/media/espace/orieux/results/';
for itask = 1:Ntasks
    
    name = [placemount,expname,'/gpac_',num2str(itask)];
    load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo');
    
    dist1(itask) = sum(sum(abs(sky(:,:,2) - xchap(:,:,2))))/sum(sum(abs(sky(:,:,2)))); 
    regul(itask) = hypers(2,2,1);
    
end

sorted1 = [regul',dist1'];
sorted1 = sortrows(sorted1,1);

lesRegs = logspace(7,15,144);
theBest1 = find(lesRegs == sorted1(find(sorted1(:,2) == min(sorted1(:,2))),1));

valBest1 = lesRegs(theBest1);

clf
semilogx(sorted1(:,1),sorted1(:,2))
hold on
grid on
ylim([0 0.31])

plot(sorted1(theBest1,1), sorted1(theBest1,2),'.', 'MarkerSize', 20);
line([sorted1(theBest1,1) sorted1(theBest1,1)],[0 sorted1(theBest1,2)])
xlim([min(sorted1(:,1)) max(sorted1(:,1))])

print('-depsc',[figplace,'/bestHyper']);

%%

ligne = 200;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

%%
%% DSP

% Window Gaussienne
sigma = 30;

NalphaW = boundAlpha(2) - boundAlpha(1);
NbetaW = boundBeta(2) - boundBeta(1);

[ALPHA BETA] = ndgrid([-(NalphaW)/2:(NalphaW)/2],[-(NbetaW)/2:(NbetaW)/2]);
window = exp(-1/2 * (ALPHA.^2 + BETA.^2)/sigma^2);

%

noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1,1)./Hrond360(1)));

deltaF = 0.004;
%% rf is reduced frequency
[dspTrue rf] = estimCircularPSD(ufft2(window.*sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspNoisefree = estimCircularPSD(ufft2(window.*noisefree(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);
dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)/Hrond360(1)),deltaF);

f1 = min(find(rf/2 >= 0.001));
f2 = min(find(rf/2 >= 0.15));

%%
%%

% lesRegs = logspace(7,15,144);
% theBest = find(lesRegs == sorted1(find(sorted1(:,2) == min(sorted1(:,2))),1));

theBest = 93;

name = [placemount,expname,'/gpac_',num2str(theBest)];
load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

clf
imagesc(xchap(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)
% title(['$\gx = 1,4\times 10^{12}$'])

print('-depsc',[figplace,'/bestCirrus']);

%%
% 
% ba = boundAlpha(1):boundAlpha(2);
% bb = boundBeta(1):boundBeta(2);
% deltaF2 = 0.006;
% [psfe rf2] = estimCircularPSD(ufft2(window.*xchap(ba,bb,2))./ufft2(window.*sky(ba,bb,2)),deltaF2);
% psfc = estimCircularPSD(ufft2(window.*fillmap(ba,bb,2)/Hrond360(1))./ufft2(window.*sky(ba,bb,2)),deltaF2);
% psft = estimCircularPSD(Hrond360(:,:,1,1)/Hrond360(1),deltaF);
% 
% 
% clf
% loglog(rf2(f1:f2)/alpha_step, psft(f1:f2),'k')
% hold on
% loglog(rf2(f1:f2)/alpha_step, psfe(f1:f2),'r')
% %loglog(rf(f1:f2)/alpha_step, psfc(f1:f2))
% grid on
% legend('True PSF', 'Effective PSF')
% xlim([0.005 0.08])
% 
% break
% 
% clf
% loglog(rf/alpha_step, psft,'k')
% hold on
% loglog(rf/alpha_step, dspProposed./dspTrue,'r')
% loglog(rf/alpha_step, dspFill./dspTrue)
% grid on
% legend('True PSF', 'Effectiv PSF', 'Co-add')
% xlim([0.005 0.11])
% %ylim([1e-17 1e-7])

%%

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)
legend('True', 'Proposed', 'Co-add')
legend('Location',[0.7 0.7 0.17 0.15])

% title(['$\gx = 1,4\times 10^{12}$'])

print('-depsc',[figplace,'/sliceBestCirrus']);

%%

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
grid on
xlim([200 600])
ylim([4 7]*1e-5)
legend('True', 'Proposed', 'Co-add')
legend('Location',[0.7 0.7 0.17 0.15])

% title(['$\gx = 1,4\times 10^{12}$'])

print('-depsc',[figplace,'/sliceBestCirrusZoom']);

%%
%% DSP

dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1): ...
                                                  boundAlpha(2), ...
                                                  boundBeta(1): ...
                                                  boundBeta(2),2)),deltaF);
bestDspNormalNoisCirrus = dspProposed;

clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')
xlim([0.005 0.11])
ylim([1e-17 1e-7])
% xlim([0 rf(f2)/alpha_step])
% ylim([1e-25 1e-5])

grid on

xlabel('f [1/arcsec]')

legend('True', 'Proposed', 'Co-add', 'Conv')
legend('Location',[0.2 0.25 0.17 0.17])
% legend('boxoff')
% title(['$\gx = 1,4\times 10^{12}$'])

print('-depsc',[figplace,'/psdBestCirrus']);


%%
%%
%% Comp madmap

clf
imagesc(madmap(:,:,2)/Hrond360(1),dynamique)
axis image; axis xy; axis off
colormap gray
colorbar

print('-depsc',[figplace,'/coaddCirrus']);

%%
%% Various Noise

theBest1High = 26;
name = [placemount,'variousNoiseCirrus','/gpac_0_01_',num2str(theBest1High)];
load(name, 'xchap');

load('/media/espace/orieux/results/FillVariousNoiseCirrus/highNoise50000', 'fillmap');

dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)/Hrond360(1)),deltaF);

clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')
xlim([0.005 0.11])
ylim([1e-17 1e-7])
% xlim([0 rf(f2)/alpha_step])
% ylim([1e-25 1e-5])
grid on

xlabel('f [1/arcsec]')

legend('True', 'Proposed', 'Co-add', 'Conv')
legend('Location',[0.2 0.25 0.17 0.17])
% legend('boxoff')

print('-depsc',[figplace,'/psdHighNoiseCirrus']);

%%
theBest1Low = 36;

name = [placemount,'variousNoiseCirrus','/gpac_0_0001_',num2str(theBest1Low)];
load(name, 'xchap');

load('/media/espace/orieux/results/FillVariousNoiseCirrus/lowNoise50000', 'fillmap');

dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)/Hrond360(1)),deltaF);

clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')
xlim([0.005 0.11])
ylim([1e-17 1e-7])
% xlim([0 rf(f2)/alpha_step])
% ylim([1e-25 1e-5])
grid on

xlabel('f [1/arcsec]')
legend('True', 'Proposed', 'Co-add', 'Conv')
legend('Location',[0.2 0.25 0.17 0.17])
% legend('boxoff')

print('-depsc',[figplace,'/psdLowNoiseCirrus']);


%
% Variance

load('/media/espace/orieux/results/variousDirty/dirtymaps')
load('/media/espace/orieux/results/fillCirrus/optim50000','fillmap')

placemount = '/media/espace/orieux/results/';
expname = 'variousRegCirrusDCT';
Ntasks = 144;

name = [placemount,expname,'/gpac_',num2str(1)];
load(name, 'xchap')
dynamique = [min(min(xchap(:,:,2))) max(max(xchap(:,:,2)))];

load([placemount,expname,'/initialization'])
placemount = '/media/espace/orieux/results/';
Ntasks = 100;

cumulant1 = zeros(Nalpha, Nbeta);
cumulant2 = zeros(Nalpha, Nbeta);

cumulant1r = zeros(size(dspProposed));
cumulant2r = zeros(size(cumulant1r));

for itask = 1:Ntasks
    disp(num2str(itask))
    for iworker = 1:8
        name = [placemount,'varianceCirrusDCT2','/sample_',num2str(itask),'_',num2str(iworker)];
        load(name, 'skySample')

        im = conv2(skySample(:,:,2),fgaussian,'same');
        imr = estimCircularPSD(ufft2(window.*im(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))), deltaF);

        cumulant1 = cumulant1 + im;
        cumulant1r = cumulant1r + imr;

        cumulant2 = cumulant2 + im.^2;
        cumulant2r = cumulant2r + imr.^2;
    end
end

cumulant1 = cumulant1/(Ntasks*8);
cumulant2 = cumulant2/(Ntasks*8);

cumulant1r = cumulant1r/(Ntasks*8);
cumulant2r = cumulant2r/(Ntasks*8);

%variance = cumulant2 - cumulant1.^2;
standartDev = sqrt(cumulant2 - cumulant1.^2);

mean2(standartDev(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)))
standartDeviation = sqrt(cumulant2 - cumulant1.^2);

standartDeviationR = sqrt(cumulant2r - cumulant1r.^2);

clf
imagesc(alpha(boundAlpha(1):boundAlpha(2)),beta(boundBeta(1): boundBeta(2)), ...
        standartDeviation, [0.5*1e-6 max(max(standartDeviation))])
axis image
axis xy; axis off
colorbar
colormap(gray)

print('-depsc',[figplace,'/variance360Cirrus']);

%% DSP

clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, bestDspNormalNoisCirrus(f1:f2) + 3*standartDeviationR(f1:f2),'k--')
loglog(rf(f1:f2)/alpha_step, bestDspNormalNoisCirrus(f1:f2) - 3*standartDeviationR(f1:f2),'k--')
% xlim([0.005 0.08])
% ylim([1e-14 1e-7])
xlim([0.005 0.11])
ylim([1e-17 1e-7])

grid on
xlabel('f [1/arcsec]')

print('-depsc',[figplace,'/psdBestCirrusStd']);

break

%%
theBest = 93;
ligne = 200;

name = [placemount,expname,'/gpac_',num2str(theBest)];
load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2))
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2) + 3*standartDev(ligne,10:Nbeta), 'k--')
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2) - 3*standartDev(ligne,10:Nbeta), 'k--')
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)

% title(['$\gx = 1,4\times 10^{12}$'])

print('-depsc',[figplace,'/sliceStandartDevCirrus']);

%%

clf
plot(beta(10:Nbeta), standartDev(ligne,10:Nbeta), 'k')
grid on
xlim([beta(10) beta(end)])

% title(['$\gx = 1,4\times 10^{12}$'])

print('-depsc',[figplace,'/sliceStandartDevCirrusOnly']);

%
%% CirrusDot

load('/media/espace/orieux/results/fillCirrusDot/coadd','coaddCirrusDot')
madmap = coaddCirrusDot;
clear coaddCirrusDot;

boundAlpha = [160 420];
boundBeta = [160 420];

placemount = '/media/espace/orieux/results/';
expname = 'variousRegCirrusDotDCT';
Ntasks = 144;

load([placemount,expname,'/initialization'])

%name = [placemount,'variousRegCirrusDotDCT','/gpac_',num2str(1)];
%load(name, 'xchap')
%dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

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

theBest = 93;

name = [placemount,expname,'/gpac_',num2str(theBest)];
load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

clf
imagesc(alpha, beta, xchap(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

print('-depsc',[figplace,'/bestCirrusDot']);

%%

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
hold on
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
grid on
ylim([4 13]*1e-5)
xlim([beta(10) beta(end)])
legend('True', 'Proposed', 'Co-add')
legend('Location',[0.7 0.7 0.17 0.15])

print('-depsc',[figplace,'/sliceBestCirrusDot']);

%%

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
hold on
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
grid on
ylim([5 12]*1e-5)
xlim([200 600])
legend('True', 'Proposed', 'Co-add')
legend('Location',[0.7 0.7 0.17 0.15])

print('-depsc',[figplace,'/sliceBestCirrusDotZoom']);

%%
%% DSP

load('/media/espace/orieux/results/fillCirrusDot/optim50000','fillmap')

noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./Hrond360(1)));

deltaF = 0.004;
%% rf is reduced frequency
[dspTrue rf] = estimCircularPSD(ufft2(window.*sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
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
legend('Location',[0.2 0.25 0.17 0.17])
%legend('boxoff')

print('-depsc',[figplace,'/psdBestCirrusDot']);

% 
%% Galaxie

load('/media/espace/orieux/results/fillGalaxie/optim500000','fillmap')
load('/media/espace/orieux/results/fillGalaxie/coadd','coaddGalaxie')
madmap = coaddGalaxie;
clear coaddGalaxie;

expname = 'variousRegGalaxieDCT';
Ntasks = 144;

load([placemount,expname,'/initialization'])

%%

clf
imagesc(madmap(:,:,2)/Hrond360(1), [0 8e-5])
axis image
axis xy; axis off
colorbar
colormap(gray)

print('-depsc',[figplace,'/coaddGalaxie'])

%%

clf
imagesc(sky(:,:,2), [0 8e-5])
axis image
axis xy; axis off
colorbar
colormap(gray)

print('-depsc',[figplace,'/galaxie'])

%%
theBest1 = 94;

ligne = 365;
ligne_mm = max(find(alpha_mm < alpha(ligne)));
load('/media/espace/orieux/results/fillGalaxie/coadd','coaddGalaxie')

lesRegs = logspace(7,15,144);
theBest = find(lesRegs == sorted1(find(sorted1(:,2) == min(sorted1(:,2))),1));

name = [placemount,expname,'/gpac_',num2str(theBest)];
load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

clf
imagesc(alpha, beta, xchap(:,:,2), [0 8]*1e-5)
axis image
axis xy; axis off
colorbar
colormap(gray)

print('-depsc',[figplace,'/bestGalaxie']);

%%

clf
plot(beta, sky(ligne,:,2), 'r')
hold on
plot(beta, xchap(ligne,:,2))
hold on
plot(beta_mm, coaddGalaxie(ligne_mm,:,2)/Hrond360(1),'k')
grid on
legend('True', 'Proposed', 'Co-add')
legend('Location',[0.7 0.7 0.17 0.15])

xlim([beta(1) beta(end)])

print('-depsc',[figplace,'/sliceBestGalaxie']);

%%

clf
plot(beta, sky(ligne,:,2), 'r')
hold on
plot(beta, xchap(ligne,:,2))
hold on
plot(beta_mm, coaddGalaxie(ligne_mm,:,2)/Hrond360(1),'k')
grid on
legend('True', 'Proposed', 'Co-add')
legend('Location',[0.7 0.7 0.17 0.15])

xlim([100 350])

print('-depsc',[figplace,'/sliceBestGalaxieZoom']);

break

%%
%% DSP

noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./Hrond360(1)));

deltaF = 0.004;
%% rf is reduced frequency
[dspTrue rf] = estimCircularPSD(ufft2(window.*sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
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

xlabel('f [1/arcsec]$')
legend('True', 'Proposed', 'Co-add', 'Conv')
legend('Location',[0.2 0.25 0.17 0.17])
%legend('boxoff')

print('-depsc',[figplace,'/psdBestGalaxie']);

%

load /media/espace/orieux/results/longRunCirrus/gpac.mat

long = size(Histo,1);

clf
loglog(Histo(:,1))
grid on

print('-depsc',[figplace,'/criterion'])

clf
loglog(Histo(5:8700,2))
grid on

print('-depsc',[figplace,'/erreurF'])

clf
loglog(Histo(1:8700,3))
grid on

print('-depsc',[figplace,'/erreurX'])

