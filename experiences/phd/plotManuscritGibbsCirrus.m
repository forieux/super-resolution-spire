clear all

% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

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
%% Gibbs Cirrus
%%

placemount = '/espace/orieux/results/'

load([placemount,'/init'])

expname = 'usmseCirrus'

alphaB = [5 555];
betaB = [7 570];

boundAlpha = [160 420];
boundBeta = [160 420];

boundBetal = [160 420];

name = [placemount,'variousRegCirrusDCT/gpac_',num2str(1)];
load(name, 'xchap')
dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

figplace = '/home/orieux/tex/manuscrit/imagesRes/gibbsCirrus/';
system(['mkdir -p ',figplace,'/eps']);

%%

load([placemount,expname,'/EAP'])
sol = conv2(xEap(:,:,2),fgaussian,'same');

true = conv2(sky(:,:,2),fgaussian,'same');

load('/espace/orieux/results/variousDirty/dirtymaps')
%load('/espace/orieux/results/fillCirrus/coadd','coaddCirrus')
load('/espace/orieux/results/fillCirrus/optim50000','fillmap')

coadded = conv2(fillmap(:,:,2),fgaussian,'same');

clf
imagesc(alpha(alphaB(1):alphaB(2)), beta(betaB(1):betaB(2)), ...
        sol(alphaB(1):alphaB(2),betaB(1):betaB(2)), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)
ylabel('$\alpha$')
xlabel('$\beta$')

l1CirrusEAP = sum2(abs(sky(:,:,2) - xEap(:,:,2)))/sum2(abs(sky(:,:,2))) 
l2CirrusEAP = sum2((sky(:,:,2) - xEap(:,:,2)).^2)/sum2((sky(:,:,2).^2)) 

name = 'eapCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

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


legend('vrai', 'mcr', 'co-add')
legend('Location',[0.6 0.6 0.32 0.25])
%legend('boxoff')

grid on
xlabel('$\beta$')

name = 'sliceCirrusMoy';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

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

clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')

xlim([0.005 0.11])
ylim([1e-17 1e-7])

grid on

xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
ylabel('$|S|^2$')
legend('vrai', 'mcr', 'co-add', 'conv')
legend('Location',[0.2 0.2 0.32 0.25])
%µlegend('boxoff')

name = 'psdEapCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

load /espace/orieux/results/usmseCirrus/samples/sample_600

clf
imagesc(alpha(alphaB(1):alphaB(2)), beta(betaB(1):betaB(2)), ...
        skySample(alphaB(1):alphaB(2),betaB(1):betaB(2),2), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)
ylabel('$\alpha$')
xlabel('$\beta$')

name = 'skySampleCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

sample = conv2(skySample(:,:,2),fgaussian,'same');

clf
imagesc(alpha(alphaB(1):alphaB(2)), beta(betaB(1):betaB(2)), ...
        sample(alphaB(1):alphaB(2),betaB(1):betaB(2)), dynamique)
axis image
axis xy; axis off
colorbar
colormap(gray)
ylabel('$\alpha$')
xlabel('$\beta$')

name = 'skySampleCirrusConv';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

iter = [1:length(gbChain(2,:))];
currentMean = cumsum(gbChain(2,:));

clf
loglog(gbChain(2,:))
hold on
loglog(gbChain(2,1),'x')
xlim([iter(1) iter(end)])

grid on

name = 'chaineGammaB';
laprint(1,[figplace,'/',name], 'width',5.5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

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

name = 'sousChaineGammaC';
laprint(1,[figplace,'/',name], 'width',5.5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
hist(gbChain(2,400:end),50)
xlim([0.97 1.03]*1e6)

name = 'histGammaB';
laprint(1,[figplace,'/',name], 'width',5.5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

meanGb = mean(gbChain(2,400:end))

stdGb = sqrt(mean(gbChain(2,400:end).^2) - meanGb^2)

%%
%%

iter = [1:length(squeeze(gxChain(1,2,:)))]';
currentMean = squeeze(cumsum(squeeze(gxChain(1,2,:))))./iter;

clf
loglog(squeeze(gxChain(1,2,:)))
hold on
loglog(gxChain(1,2,1),'x')
xlim([iter(1) iter(end)])
grid on

name = 'chaineGammaX';
laprint(1,[figplace,'/',name], 'width',5.5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

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

name = 'sousChaineGammaX';
laprint(1,[figplace,'/',name], 'width',5.5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

meanGx = mean(squeeze(gxChain(1,2,400:end)))
stdGx = sqrt(mean(gxChain(1,2,400:end).^2) - meanGx^2)

%%

clf
hist(gxChain(1,2,400:end),50)
xlim([2.6 2.9]*1e11)

name = 'histGammaX';
laprint(1,[figplace,'/',name], 'width',5.5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% Variance

% cumulant1 = zeros(size(sample));
% cumulant2 = zeros(size(sample));
% hit = 0;

% for itask = 400:932
    
%     disp(num2str(itask));
    
%     name = [placemount,expname,'/samples/sample_',num2str(itask)];
%     try
%         load(name, 'skySample')
        
%         cumulant1 = cumulant1 + conv2(skySample(:,:,2),fgaussian,'same');
%         cumulant2 = cumulant2 + (conv2(skySample(:,:,2),fgaussian, 'same')).^2;
        
%         hit = hit+1;
%     end
% end

% cumulant1 = cumulant1/hit;
% cumulant2 = cumulant2/hit;
% standartDev = sqrt(cumulant2 - cumulant1.^2);

% mean2(standartDev(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)))

% clf
% imagesc(alpha,beta,sqrt(cumulant2 - cumulant1.^2))
% axis image
% axis xy; axis off
% colorbar
% colormap(gray)
% ylabel('$\beta$')
% xlabel('$\alpha$')

% name = 'varianceCirrus';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

% clf
% plot(beta(10:end), sol(ligne,10:end))
% hold on
% plot(beta(10:end), sol(ligne,10:end) + 3*standartDev(ligne,10:end), 'k')
% plot(beta(10:end), sol(ligne,10:end) - 3*standartDev(ligne,10:end), 'k')
% grid on
% xlim([beta(10) beta(end)])

% xlabel('$\beta$')

% name = 'sliceStandartDevGibbsCirrus';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

close all
