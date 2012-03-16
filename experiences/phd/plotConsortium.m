clear all

% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

figplace = '/home/orieux/tex/planches/workshop/';
system(['mkdir -p ',figplace]);

set(0,'defaulttextinterpreter','tex')

set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxeslinewidth',1);
set(0,'defaultlinelinewidth',2);

SizeMarker = 2;
SizeLabel = 30;

numfig = 1;
figure(numfig)

%%
%% Cirrus
%%

load('/espace/orieux/consortium/dirtymaps')
load('/espace/orieux/consortium/optim50000','fillmap')
load('/espace/orieux/consortium/initialization')

placemount = '/espace/orieux/'
expname = 'consortium';
Ntasks = 144;

boundAlpha = [160 420];
boundBeta = [160 420];

boundAlpha2 = [465 950];
boundBeta2 = [-85 465];
 
boundBetal = [160 420];

name = [placemount,expname,'/gpac_93'];
load(name, 'xchap')
dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

load([placemount,expname,'/initialization'])

%%

clf
imagesc(alpha, beta, sky(:,:,2), dynamique)
axis equal
axis xy; axis off
colorbar
colormap(gray)
% ylim(boundBeta2)
% xlim(boundAlpha2)

name = 'trueCirrus';
print(1,'-depsc',[figplace,'/',name])

ligne = 200;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

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

%

theBest = 93;

clf
imagesc(alpha, beta, xchap(:,:,2), dynamique)
axis equal
axis xy; axis off
colorbar
colormap(gray)
% ylim(boundBeta2)
% xlim(boundAlpha2)

name = 'bestCirrus';
print(1,'-depsc',[figplace,'/',name])

%%
%% DSP

noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./Hrond360(1)));

deltaF = 0.004;
[dspTrue rf] = estimCircularPSD(ufft2(window.*sky(boundAlpha(1): ...
                                                  boundAlpha(2), boundBeta(1): ...
                                                  boundBeta(2),2)),deltaF);
dspNoisefree = estimCircularPSD(ufft2(window.*noisefree(boundAlpha(1): ...
                                                  boundAlpha(2), boundBeta(1): ...
                                                  boundBeta(2))),deltaF);
dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1): boundAlpha(2), ...
                                                  boundBeta(1): boundBeta(2),2)/ ...
                                  Hrond360(1)),deltaF);

f1 = min(find(rf/2 >= 0.001));
f2 = min(find(rf/2 >= 0.15));


dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1): ...
                                                  boundAlpha(2), boundBeta(1): ...
                                                  boundBeta(2),2)),deltaF);

clf
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspProposed(f1:f2))
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspFill(f1:f2),'k')
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspNoisefree(f1:f2),'--k')
xlim([0.005 0.11])
ylim([1e-17 1e-7])
% xlim([0 rf(f2)/alpha_step])
% ylim([1e-25 1e-5])

grid on

xlabel('arcsec^{-1}')
ylabel('|S|^2')
legend('true', 'proposed', 'co-add', 'conv')
%legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')
% title(['\gx = 1,4\times 10^{12}'])

name = 'psdBestCirrus';
print(1,'-depsc',[figplace,'/',name])

%%

deltaF = 0.001;
Hrond = (abs(Hrond360(:,:,1,1)).^2 + ((1.4*10^12)/10^6)*abs(circDalpha + ...
                                                  circDbeta).^2).^(-1).* ...
        abs(Hrond360(:,:,1,1)).^2;
[HmeanCorrected rf] = estimCircularPSD(Hrond,deltaF);
Hmean = estimCircularPSD(Hrond360(:,:,1,1),deltaF);

clf
loglog(rf*sqrt(1/alpha_step^2+1/beta_step^2),HmeanCorrected./max(HmeanCorrected))
hold on
loglog(rf*sqrt(1/alpha_step^2+1/beta_step^2),Hmean./max(Hmean),'--')

legend('proposed effective MTF', 'co-add effective MTF','Location', 'SW')
    
grid on
xlim([0.003 0.11]) 
ylim([1e-8 1e1])

xlabel('arcsec^{-1}','interpreter','tex')

name = 'effectivePSF';
print(1,'-depsc',[figplace,'/',name])


clf
plot(rf*sqrt(1/alpha_step^2+1/beta_step^2),HmeanCorrected./max(HmeanCorrected))
hold on
plot(rf*sqrt(1/alpha_step^2+1/beta_step^2),Hmean./max(Hmean),'--')

% falpha = linspace(-0.5,0.5,Nalpha);
% fbeta = linspace(-0.5,0.5,Nbeta);

% [FALPHAS FBETAS] = ndgrid(falpha,fbeta);
% FALPHAS = circshift(ifftshift(FALPHAS), [floor(Nsalpha/2) 1]);
% FBETAS = circshift(ifftshift(FBETAS), [1 floor(Nsbeta/2)]);

% H = circshift(Hrond360(:,:,1,1), [floor(Nsalpha/2) floor(Nsbeta/2)]);

% plot(FALPHAS(1:Nsalpha,1)/alpha_step, abs(H(1:Nsalpha, floor(Nsbeta/2)))/max2(abs(H(1:Nsalpha, floor(Nsbeta/2)))),'r--')

legend('proposed effective MTF', 'co-add effective MTF')
    
grid on
xlim([0 0.07]) 
%ylim([1e-8 1e1])

xlabel('arcsec^{-1}','interpreter','tex')

name = 'effectivePSFlin';
print(1,'-depsc',[figplace,'/',name])


%%
%% Comp madmap

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)
legend('true', 'co-add', 'proposed')
%legend('Location',[0.6 0.6 0.32 0.25])

xlabel('arcsec')

name = 'sliceCompCoaddDeconv';
print(1,'-depsc',[figplace,'/',name])

%%

clf
imagesc(alpha_mm, beta_mm, madmap(:,:,2)/Hrond360(1), dynamique)
axis equal
axis xy; axis off
colorbar
colormap(gray)
% ylim(boundBeta2)
% xlim(boundAlpha2)


name = 'coadd360cirrus';
print(1,'-depsc',[figplace,'/',name])

%%

% load('/espace/orieux/results/variousDirty/dirtymaps')
% load('/espace/orieux/results/fillCirrus/optim50000','fillmap')

% placemount = '/espace/orieux/results/'
% expname = 'variousRegCirrusDCT';
% Ntasks = 144;


% name = [placemount,expname,'/gpac_',num2str(1)];
% load(name, 'xchap')
% dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

% load([placemount,expname,'/initialization'])

% Ntasks = 100;

% cumulant1 = zeros(Nalpha, Nbeta);
% cumulant2 = zeros(Nalpha, Nbeta);

% cumulant1r = zeros(NalphaW+1, NbetaW+1);
% cumulant2r = zeros(NalphaW+1, NbetaW+1);

% for itask = 1:Ntasks
%     disp(num2str(itask))
%     for iworker = 1:8
%         name = [placemount,'varianceCirrusDCT2','/sample_',num2str(itask),'_',num2str(iworker)];
%         load(name, 'skySample')
        
%         im = conv2(skySample(:,:,2),fgaussian,'same');
%         imr = ufft2(window.*im(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)));

%         cumulant1 = cumulant1 + im;
%         cumulant1r = cumulant1r + imr;

%         cumulant2 = cumulant2 + im.^2;
%         cumulant2r = cumulant2r + imr.^2;
%     end
% end

% cumulant1 = cumulant1/(Ntasks*8);
% cumulant2 = cumulant2/(Ntasks*8);

% cumulant1r = cumulant1r/(Ntasks*8);
% cumulant2r = cumulant2r/(Ntasks*8);

% %variance = cumulant2 - cumulant1.^2;
% standartDev = sqrt(cumulant2 - cumulant1.^2);

% mean2(standartDev(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)))
% standartDeviation = sqrt(cumulant2 - cumulant1.^2);

% standartDeviationR = sqrt(cumulant2r - cumulant1r.^2);

% clf
% imagesc(alpha(boundAlpha(1):boundAlpha(2)),beta(boundBeta(1): boundBeta(2)), ...
%         standartDeviation, [0.5*1e-6 max2(standartDeviation)])
% axis equal
% axis xy; axis off
% colorbar
% colormap(gray)

% name = 'variance360Cirrus';
% print(1,'-depsc',[figplace,'/',name])

% %%
% theBest = 93;
% ligne = 200;

% name = [placemount,expname,'/gpac_',num2str(theBest)];
% load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
% hold on
% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2) + 3*standartDev(ligne,10:Nbeta), 'k')
% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2) - 3*standartDev(ligne,10:Nbeta), 'k')
% grid on
% xlim([beta(10) beta(end)])
% ylim(dynamique)

% xlabel('arcsec')
% % title(['\gx = 1,4\times 10^{12}'])

% name = 'sliceStandartDevCirrus';
% print(1,'-depsc',[figplace,'/',name])

%

% noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./Hrond360(1)));

% load('/espace/orieux/results/fillCirrus/optim50000','fillmap')

% deltaF = 0.004;
% %% rf is reduced frequency
% [dspTrue rf] = estimCircularPSD(ufft2(window.*sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
% dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
% dspSTD = estimCircularPSD(standartDeviationR,deltaF);
% dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)/Hrond360(1)),deltaF);


% f1 = min(find(rf/2 >= 0.001));
% f2 = min(find(rf/2 >= 0.15));

% clf
% % RF = F/FE = F*Te

% loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspTrue(f1:f2),'r')
% hold on
% loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspProposed(f1:f2))
% loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspFill(f1:f2),'k')
% loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspProposed(f1:f2)+3*dspSTD(f1:f2),'--')
% legend('true', 'proposed', 'coadd', 'std')
% %legend('Location',[0.2 0.2 0.32 0.25])
% loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspProposed(f1:f2)-3*dspSTD(f1:f2),'--')

% xlim([0.005 0.11])
% ylim([1e-17 1e-7])
% % xlim([0.005 0.11])
% % ylim([1e-17 1e-7])
% grid on

% xlabel('{arcsec}^{-1}')
% ylabel('|S|^2')
% % legend('boxoff')

% % title(['\gx = 10^7'])

% name = 'psdStandartDevCirrus';
% print(1,'-depsc',[figplace,'/',name])

%%
%% Galaxie
%%


load('/espace/orieux/consortium/optim500000_gal','fillmap')
load('/espace/orieux/consortium/initialization_gal')

ligne = 365;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

clf
imagesc(alpha, beta, sky(:,:,2), [0 8e-5])
axis equal
axis xy; axis off
colorbar
colormap(gray)
% ylim(boundBeta2)
% xlim(boundAlpha2)

name = 'trueGal';
print(1,'-depsc',[figplace,'/',name])

%%
%% Comp madmap

load('/espace/orieux/consortium/coadd_gal','coaddGalaxie')

theBest = 94;
load('/espace/orieux/consortium/gpac_94_gal')


clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta_mm, coaddGalaxie(ligne_mm,:,2)/Hrond360(1),'k')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
grid on
% xlim([beta(10) beta(end)])
% ylim(dynamique)
legend('true', 'co-add', 'proposed')
%legend('Location',[0.6 0.6 0.32 0.25])

xlabel('arcsec')

name = 'sliceCompCoaddDeconvGal';
print(1,'-depsc',[figplace,'/',name])

%%


clf
imagesc(alpha_mm, beta_mm, coaddGalaxie(:,:,2)/Hrond360(1), [0 8e-5])
axis equal
axis xy; axis off
colorbar
colormap(gray)
% ylim(boundBeta2)
% xlim(boundAlpha2)

name = 'coadd360gal';
print(1,'-depsc',[figplace,'/',name])

%%

ligne = 365;
ligne_mm = max(find(alpha_mm < alpha(ligne)));



clf
imagesc(alpha, beta, xchap(:,:,2), [0 8]*1e-5)
axis image
axis xy; axis off
colorbar
colormap(gray)
ylabel('\alpha')
xlabel('\beta')

name = 'bestGalaxie';
print(1,'-depsc',[figplace,'/',name])

%%

clf
plot(beta, sky(ligne,:,2), 'r')
hold on
plot(beta, xchap(ligne,:,2))
hold on
plot(beta_mm, coaddGalaxie(ligne_mm,:,2)/Hrond360(1),'k')
grid on
legend('true', 'proposed', 'co-add')
%legend('Location',[0.6 0.6 0.32 0.25])

xlim([beta(1) beta(end)])

xlabel('arcsec')

name = 'sliceBestGalaxie';
print(1,'-depsc',[figplace,'/',name])

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
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspProposed(f1:f2))
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspFill(f1:f2),'k')
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspNoisefree(f1:f2),'--k')

xlim([0.005 0.11])
ylim([1e-17 1e-6])

grid on

xlabel('arcsec^{-1}')
ylabel('|S|^2')
legend('true', 'proposed', 'co-add', 'conv')
%legend('Location',[0.2 0.2 0.32 0.25])
%legend('boxoff')

name = 'psdBestGalaxie';
print(1,'-depsc',[figplace,'/',name])

%%
%%
%% CirrusDot
%%

ligne = 200;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

load('/espace/orieux/consortium/optim50000_cd','fillmap')

boundAlpha = [160 420];
boundBeta = [160 420];

Ntasks = 144;

load('/espace/orieux/consortium/initialization_cd')
load('/espace/orieux/consortium/gpac_93_cd')

load('/espace/orieux/consortium/coadd_cd')

%%

clf
imagesc(alpha, beta, sky(:,:,2), dynamique)
axis equal
axis xy; axis off
colorbar
colormap(gray)
% ylim(boundBeta2)
% xlim(boundAlpha2)

name = 'trueCirrusDot';
print(1,'-depsc',[figplace,'/',name])

%%
%% Comp madmap

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta_mm, coaddCirrusDot(ligne_mm,:,2)/Hrond360(1),'k')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)
legend('true', 'co-add', 'proposed')
%legend('Location',[0.6 0.6 0.32 0.25])

xlabel('arcsec')

name = 'sliceCompCoaddDeconvCirrusDot';
print(1,'-depsc',[figplace,'/',name])

%%

clf
imagesc(alpha_mm, beta_mm, coaddCirrusDot(:,:,2)/Hrond360(1), dynamique)
axis equal
axis xy; axis off
colorbar
colormap(gray)
% ylim(boundBeta2)
% xlim(boundAlpha2)

name = 'coadd360CirrusDot';
print(1,'-depsc',[figplace,'/',name])

%%

dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

ligne = 340;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

theBest = 93;

clf
imagesc(alpha, beta, xchap(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

name = 'bestCirrusDot';
print(1,'-depsc',[figplace,'/',name])

%%

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
hold on
plot(beta_mm, coaddCirrusDot(ligne_mm,:,2)/Hrond360(1),'k')
grid on
ylim([4 13]*1e-5)
xlim([beta(10) beta(end)])
legend('true', 'proposed', 'co-add')
%legend('Location',[0.6 0.6 0.32 0.25])

xlabel('arcsec')

name = 'sliceBestCirrusDot';
print(1,'-depsc',[figplace,'/',name])

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
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspProposed(f1:f2))
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspFill(f1:f2),'k')
loglog(rf(f1:f2)*sqrt(1/alpha_step^2+1/beta_step^2), dspNoisefree(f1:f2),'--k')

xlim([0.005 0.11])
ylim([1e-17 1e-7])

grid on

xlabel('arcsec^{-1}')
ylabel('|S|^2')
legend('true', 'proposed', 'co-add', 'conv')
%legend('Location',[0.2 0.2 0.32 0.25])
%legend('boxoff')

name = 'psdBestCirrusDot';
print(1,'-depsc',[figplace,'/',name])

close all
