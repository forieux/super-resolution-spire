clear all

% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

figplace = '/home/orieux/tex/manuscrit/imagesRes/cirrus/';
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
expname = 'variousRegCirrusDCT';
Ntasks = 144;

boundAlpha = [160 420];
boundBeta = [160 420];

boundBetal = [160 420];


name = [placemount,expname,'/gpac_',num2str(1)];
load(name, 'xchap')
dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

load([placemount,expname,'/initialization'])

set(0,'defaulttextinterpreter','none')

%%

clf
imagesc(sky(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

name = 'true';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%% L2 et L1

for itask = 1:Ntasks
    
    disp(num2str(itask));
    
    name = [placemount,expname,'/gpac_',num2str(itask)];
    load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')
    
    dist1(itask) = sum2(abs(sky(:,:,2) - xchap(:,:,2)))/sum2(abs(sky(:,:,2))); 
    dist2(itask) = sum2((sky(:,:,2) - xchap(:,:,2)).^2)/sum2((sky(:,:,2).^2)); 
    regul(itask) = hypers(2,2,1);
    
end

sorted1 = [regul',dist1'];
sorted1 = sortrows(sorted1,1);
sorted2 = [regul',dist2'];
sorted2 = sortrows(sorted2,1);

lesRegs = logspace(7,15,144);
theBest1 = find(lesRegs == sorted1(find(sorted1(:,2) == min(sorted1(:,2))),1));
theBest2 = find(lesRegs == sorted2(find(sorted2(:,2) == min(sorted2(:,2))),1));

valBest1 = lesRegs(theBest1)
valBest2 = lesRegs(theBest2)

sorted1(theBest1,2)
sorted2(theBest2,2)

clf
semilogx(sorted1(:,1),sorted1(:,2))
hold on
semilogx(sorted2(:,1),sorted2(:,2),'r')
grid on
ylim([0 0.31])
legend('$\ell_1$','$\ell_2$', 'Location', [0.5 0.7 0.3 0.2])
legend('boxoff')

plot(sorted1(theBest1,1), sorted1(theBest1,2),'.', 'MarkerSize', 10);
plot(sorted2(theBest2,1), sorted2(theBest2,2),'.r', 'MarkerSize', 10);

line([sorted1(theBest1,1) sorted1(theBest1,1)],[0 sorted1(theBest1,2)])
h = line([sorted2(theBest2,1) sorted2(theBest2,1)],[0 sorted2(theBest2,2)]);
set(h,'Color',[1 0 0]);

xlim([min(sorted1(:,1)) max(sorted1(:,1))])

xlabel('$\gx~(\gB = 10^6)$')

name = 'bestHyper';
laprint(1,[figplace,'/',name], 'width',7)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

break

%%
%%

name = [placemount,expname,'/gpac_',num2str(1)];
load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

clf
imagesc(xchap(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)
% title(['$\gx = 10^7$']);%,num2str(hypers(2,1,1),'%0.1e'),'$'])

name = 'lowRegCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

ligne = 200;
ligne_mm = max(find(alpha_mm < alpha(ligne)));

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)
legend('vrai', 'mcr', 'co-add')
legend('Location',[0.6 0.6 0.32 0.25])

xlabel('$\beta $')
% title(['$\gx = 10^7$'])

name = 'slicelowRegCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

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
dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
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
xlim([0.005 0.11])
ylim([1e-17 1e-7])
% xlim([0.005 0.11])
% ylim([1e-17 1e-7])
grid on

xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
ylabel('$|S|^2$')
legend('vrai', 'mcr', 'co-add', 'conv')
legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')

% title(['$\gx = 10^7$'])

name = 'psdLowRegCirrus';
laprint(1,[figplace,'/',name], 'width',8, 'asonscreen', 'on')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

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

name = 'bestCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)
legend('vrai', 'mcr', 'co-add')
legend('Location',[0.6 0.6 0.32 0.25])

xlabel('$\beta$')
% title(['$\gx = 1,4\times 10^{12}$'])

name = 'sliceBestCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% DSP

dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1): ...
                                                  boundAlpha(2), ...
                                                  boundBeta(1): ...
                                                  boundBeta(2),2)),deltaF);

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

xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
ylabel('$|S|^2$')
legend('vrai', 'mcr', 'co-add', 'conv')
legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')
% title(['$\gx = 1,4\times 10^{12}$'])

name = 'psdBestCirrus';
laprint(1,[figplace,'/',name], 'width',8, 'asonscreen', 'on')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])


% clf
% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
% hold on

% name = [placemount,expname,'/gpac_',num2str(theBest)];
% load(name, 'xchap')

% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2),'k')
% grid on
% xlabel('$\beta$')

% name = 'sliceLongRun';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

name = [placemount,expname,'/gpac_',num2str(Ntasks)];
load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

clf
imagesc(xchap(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)
% title(['$\gx = 10^{15}$'])

name = 'tomuchRegCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)
legend('vrai', 'mcr', 'co-add')
legend('Location',[0.6 0.6 0.32 0.25])

xlabel('$\beta$')
% title(['$\gx = 10^{15}$'])

name = 'sliceTomuchRegCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% DSP

dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);

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

xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
ylabel('$|S|^2$')
legend('vrai', 'mcr', 'co-add', 'conv')
legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')
% title(['$\gx = 10^{15}$'])

name = 'psdTomuchRegCirrus';
laprint(1,[figplace,'/',name], 'width',8, 'asonscreen', 'on')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
load /espace/orieux/results/longRunCirrus/gpac.mat

long = size(Histo,1);

clf
loglog(Histo(:,1))
grid on
title('$J\left(\xb^{(q)}\right)$')

name = 'criterion';
laprint(1,[figplace,'/',name], 'width',5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

clf
loglog(Histo(5:8700,2))
grid on
title('$\left\|J\left(\xb^{(q)}\right) - J\left(\xb^{(q-1)}\right)\right\|^2$')

name = 'erreurF';
laprint(1,[figplace,'/',name], 'width',5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

clf
loglog(Histo(1:8700,3))
grid on
title('$\left\|\xb^{(q)} - \xb^{(q-1)}\right\|^2$')

name = 'erreurX';
laprint(1,[figplace,'/',name], 'width',5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

clf
imagesc(xchap(:,:,2), dynamique)
axis image
axis xy; axis off
axis off
colorbar
colormap(gray)

name = 'longRun';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%% Comp madmap

name = [placemount,expname,'/gpac_',num2str(theBest)];
load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2), 'r')
hold on
plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)
legend('vrai', 'mcr', 'co-add')
legend('Location',[0.6 0.6 0.32 0.25])

xlabel('$\beta$')

name = 'sliceCompCoaddDeconv';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
imagesc(alpha_mm, beta_mm, madmap(:,:,2)/Hrond360(1), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)

name = 'coadd360cirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%% Various Noise

for itask = 1:40
    
    disp(num2str(itask));
    
    name = [placemount,'variousNoiseCirrus','/gpac_0_01_',num2str(itask)];
    load(name, 'xchap', 'hypers');%, 'SortieOPS', 'Histo')
    

    dist1High(itask) = sum2(abs(sky(boundAlpha(1):boundAlpha(2), ...
                                    boundBeta(1):boundBeta(2),2) - ...
                                xchap(boundAlpha(1):boundAlpha(2), ...
                                      boundBeta(1):boundBeta(2),2)))/ ...
        sum2(abs(sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)));         
    
    dist2High(itask) = sum2((sky(boundAlpha(1):boundAlpha(2), ...
                                 boundBeta(1):boundBeta(2),2) - ...
                             xchap(boundAlpha(1):boundAlpha(2), ...
                                   boundBeta(1):boundBeta(2),2)).^2)/ ...
        sum2((sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2).^2)); 
    
    regulHigh(itask) = hypers(2,1,1);
end

sorted1High = [regulHigh',dist1High'];
sorted1High = sortrows(sorted1High,1);
sorted2High = [regulHigh',dist2High'];
sorted2High = sortrows(sorted2High,1);

lesRegs = logspace(10,13,40);
theBest1High = find(lesRegs == sorted1High(find(sorted1High(:,2) == min(sorted1High(:,2))),1))
theBest2High = find(lesRegs == sorted2High(find(sorted2High(:,2) == min(sorted2High(:,2))),1))

valBest1High = lesRegs(theBest1High)
valBest2High = lesRegs(theBest2High)

plot(sorted1High(:,1), sorted1High(:,2))

name = [placemount,'variousNoiseCirrus','/gpac_0_01_',num2str(theBest1High)];
load(name, 'xchap');

clf
imagesc(xchap(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)
% title('$\gB = 10^{4}$ et $\gx = 8,3 \times 10^{11}$')

name = 'highNoiseCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

load('/espace/orieux/results/FillVariousNoiseCirrus/highNoise50000', 'fillmap');

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2),'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
plot(beta(10:Nbeta), fillmap(ligne,10:Nbeta,2)/Hrond360(1),'k')
xlabel('$\beta$')
legend('vrai', 'mcr', 'co-add')
legend('Location',[0.6 0.6 0.32 0.25])
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)

name = 'slicehighNoiseCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

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

xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
ylabel('$|S|^2$')
legend('vrai', 'mcr', 'co-add', 'conv')
legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')

name = 'psdHighNoiseCirrus';
laprint(1,[figplace,'/',name], 'width',8, 'asonscreen', 'on')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

for itask = 1:40
    
    disp(num2str(itask));
    
    name = [placemount,'variousNoiseCirrus','/gpac_0_0001_',num2str(itask)];
    load(name, 'xchap', 'hypers');%, 'SortieOPS', 'Histo')
    
    dist1Low(itask) = sum2(abs(sky(boundAlpha(1):boundAlpha(2), ...
                                    boundBeta(1):boundBeta(2),2) - ...
                                xchap(boundAlpha(1):boundAlpha(2), ...
                                      boundBeta(1):boundBeta(2),2)))/ ...
        sum2(abs(sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)));         
    
    dist2Low(itask) = sum2((sky(boundAlpha(1):boundAlpha(2), ...
                                 boundBeta(1):boundBeta(2),2) - ...
                             xchap(boundAlpha(1):boundAlpha(2), ...
                                   boundBeta(1):boundBeta(2),2)).^2)/ ...
        sum2((sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2).^2)); 
 
    regulLow(itask) = hypers(2,2,1);
end

sorted1Low = [regulLow',dist1Low'];
sorted1Low = sortrows(sorted1Low,1);
sorted2Low = [regulLow',dist2Low'];
sorted2Low = sortrows(sorted2Low,1);

lesRegs = logspace(10,13,40);
theBest1Low = find(lesRegs == sorted1Low(find(sorted1Low(:,2) == min(sorted1Low(:,2))),1))
theBest2Low = find(lesRegs == sorted2Low(find(sorted2Low(:,2) == min(sorted2Low(:,2))),1))

valBest1Low = lesRegs(theBest1Low)
valBest2Low = lesRegs(theBest2Low)

name = [placemount,'variousNoiseCirrus','/gpac_0_0001_',num2str(theBest1Low)];
load(name, 'xchap');

clf
imagesc(xchap(:,:,2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)
% title('$\gB = 10^{8}$ et  $\gx = 4,9 \times 10^{12}$')

name = 'lowNoiseCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

load('/espace/orieux/results/FillVariousNoiseCirrus/lowNoise50000', 'fillmap');

clf
plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2),'r')
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
plot(beta(10:Nbeta), fillmap(ligne,10:Nbeta,2)/Hrond360(1),'k')
xlabel('$\beta$')
legend('vrai', 'mcr', 'co-add')
legend('Location',[0.6 0.6 0.32 0.25])
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)

name = 'slicelowNoiseCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

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

xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
ylabel('$|S|^2$')
legend('vrai', 'mcr', 'co-add', 'conv')
legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')

name = 'psdLowNoiseCirrus';
laprint(1,[figplace,'/',name], 'width',8, 'asonscreen', 'on')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%%

% clf
% plot(beta(10:Nbeta), fillmap(ligne,10:Nbeta,2)/Hrond360(1), 'r')
% hold on
% plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
% grid on
% xlim([beta(10) beta(end)])

% xlabel('$\beta$')

% name = 'sliceCompCoaddFill';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
% %% Three band

% load /espace/orieux/results/cirrusThreeBandIndep/gpac

% clf
% imagesc(xchap(:,:,1), dynamique)
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')
% title('$250~\mu\textrm{m}$')

% name = 'cirrus250';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

% clf
% plot(beta(10:Nbeta), sky(ligne,10:Nbeta,1),'r')
% hold on
% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,1))
% xlabel('$\beta$')
% title('$250~\mu\textrm{m}$')
% grid on

% name = 'slicecirrus250';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

% clf
% imagesc(xchap(:,:,2), dynamique)
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')
% title('$360~\mu\textrm{m}$')

% name = 'cirrus360';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

% clf
% plot(beta(10:Nbeta), sky(ligne,10:Nbeta,2),'r')
% hold on
% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
% xlabel('$\beta$')
% title('$360~\mu\textrm{m}$')
% grid on

% name = 'slicecirrus360';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

% clf
% imagesc(xchap(:,:,3), dynamique)
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')
% title('$520~\mu\textrm{m}$')

% name = 'cirrus520';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

% clf
% plot(beta(10:Nbeta), sky(ligne,10:Nbeta,3),'r')
% hold on
% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,3))
% xlabel('$\beta$')
% title('$520~\mu\textrm{m}$')
% grid on

% name = 'slicecirrus520';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

% clf
% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,1))
% hold on
% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2),'k')
% plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,3),'g')
% xlabel('$\beta$')
% grid on

% name = 'slicecirrusThreeBand';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

% clf
% imagesc(mapmean(:,:,1)/Hrond250(1), dynamique)
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)

% name = 'dirtyMean250';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% clf
% imagesc(mapmean(:,:,2)/Hrond360(1), dynamique)
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)

% name = 'dirtyMean360';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% clf
% imagesc(mapmean(:,:,3)/Hrond520(1), dynamique)
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)

% name = 'dirtyMean520';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

% clf
% imagesc(fillmap(:,:,1)/Hrond250(1), dynamique)
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)

% name = 'fillMap250';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% clf
% imagesc(fillmap(:,:,2)/Hrond360(1), dynamique)
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)

% name = 'fillMap360';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% clf
% imagesc(fillmap(:,:,3)/Hrond520(1), dynamique)
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)

% name = 'fillMap520';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])
