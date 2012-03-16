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

boundAlpha = [160 420];
boundBeta = [160 420];

boundBetal = [160 420];



%%
%% DSP

% Window

%% Gaussienne
sigma = 30;

NalphaW = boundAlpha(2) - boundAlpha(1);
NbetaW = boundBeta(2) - boundBeta(1);

[ALPHA BETA] = ndgrid([-(NalphaW)/2:(NalphaW)/2],[-(NbetaW)/2:(NbetaW)/2]);
window = exp(-1/2 * (ALPHA.^2 + BETA.^2)/sigma^2);


%%
%% Cirrus
%%

load('/espace/orieux/results/variousDirty/dirtymaps')
load('/espace/orieux/results/fillCirrus/optim50000','fillmap')

placemount = '/espace/orieux/results/'
expname = 'variousRegCirrusDCT';
Ntasks = 144;


name = [placemount,expname,'/gpac_',num2str(1)];
load(name, 'xchap')
dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

load([placemount,expname,'/initialization'])

set(0,'defaulttextinterpreter','none')

Ntasks = 100;

cumulant1 = zeros(Nalpha, Nbeta);
cumulant2 = zeros(Nalpha, Nbeta);

cumulant1r = zeros(NalphaW+1, NbetaW+1);
cumulant2r = zeros(NalphaW+1, NbetaW+1);

for itask = 1:Ntasks
    disp(num2str(itask))
    for iworker = 1:8
        name = [placemount,'varianceCirrusDCT2','/sample_',num2str(itask),'_',num2str(iworker)];
        load(name, 'skySample')
        
        im = conv2(skySample(:,:,2),fgaussian,'same');
        imr = ufft2(window.*im(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)));

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
        standartDeviation, [0.5*1e-6 max2(standartDeviation)])
axis image
axis xy; axis off
colorbar
colormap(gray)

break

name = 'variance360Cirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
theBest = 93;
ligne = 200;

name = [placemount,expname,'/gpac_',num2str(theBest)];
load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2))
hold on
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2) + 3*standartDev(ligne,10:Nbeta), 'k')
plot(beta(10:Nbeta), xchap(ligne,10:Nbeta,2) - 3*standartDev(ligne,10:Nbeta), 'k')
grid on
xlim([beta(10) beta(end)])
ylim(dynamique)

xlabel('$\beta$')
% title(['$\gx = 1,4\times 10^{12}$'])

name = 'sliceStandartDevCirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%

noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./Hrond360(1)));

load('/espace/orieux/results/fillCirrus/optim50000','fillmap')

deltaF = 0.004;
%% rf is reduced frequency
[dspTrue rf] = estimCircularPSD(ufft2(window.*sky(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspProposed = estimCircularPSD(ufft2(window.*xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)),deltaF);
dspSTD = estimCircularPSD(standartDeviationR,deltaF);
dspFill = estimCircularPSD(ufft2(window.*fillmap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),2)/Hrond360(1)),deltaF);


f1 = min(find(rf/2 >= 0.001));
f2 = min(find(rf/2 >= 0.15));

clf
% RF = F/FE = F*Te

loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2)+3*dspSTD(f1:f2),'--')
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2)-3*dspSTD(f1:f2),'--')

xlim([0.005 0.11])
ylim([1e-17 1e-7])
% xlim([0.005 0.11])
% ylim([1e-17 1e-7])
grid on

xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
ylabel('$|S|^2$')
legend('vrai', 'mcr', 'coadd', 'std')
legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')

% title(['$\gx = 10^7$'])

