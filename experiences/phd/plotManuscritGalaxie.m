% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

figplace = '/home/orieux/tex/manuscrit/imagesRes/galaxie/';
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
%% Galaxie
%%

load('/espace/orieux/results/fillGalaxie/optim500000','fillmap')

placemount = '/espace/orieux/results/'
expname = 'variousRegGalaxieDCT';
Ntasks = 144;

load([placemount,expname,'/initialization'])

%%

% mapmean = dirtymap(data, [1 1 1], index250, index360, index520, coefs250, ...
%                    coefs360, coefs520, Nalpha, Nbeta, 'm');

% clf
% imagesc(alpha, beta, mapmean(:,:,1)/Hrond250(1))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')

% name = 'dirtyMean250Galaxie';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% clf
% imagesc(alpha, beta, mapmean(:,:,2)/Hrond360(1))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')

% name = 'dirtyMean360Galaxie';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% clf
% imagesc(alpha, beta, mapmean(:,:,3)/Hrond520(1))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')

% name = 'dirtyMean520Galaxie';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%
% %%

% clf
% imagesc(alpha, beta, fillmap(:,:,1)/Hrond250(1))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')

% name = 'fillGalaxie250';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% clf
% imagesc(alpha, beta, fillmap(:,:,2)/Hrond360(1))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')

% name = 'fillGalaxie360';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% clf
% imagesc(alpha, beta, fillmap(:,:,3)/Hrond520(1))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')

% name = 'fillGalaxie520';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

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

% clf
% semilogx(sorted1(:,1),sorted1(:,2))
% hold on
% semilogx(sorted2(:,1),sorted2(:,2),'r')
% grid on
% legend('$\ell_1$','$\ell_2$')
% plot(sorted1(theBest1,1), sorted1(theBest1,2),'.');%, 'MarkerSize', SizeMarker)
% plot(sorted2(theBest2,1), sorted2(theBest2,2),'.r');%, 'MarkerSize', SizeMarker)

% xlim([min(sorted1(:,1)) max(sorted1(:,1))])

% xlabel('$\gxb$')

% name = 'bestHyperGalaxie';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%
% %%

% name = [placemount,expname,'/gpac_',num2str(1)];
% load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

% clf
% imagesc(alpha, beta, xchap(:,:,2))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')
% title(['$\gxb = 10^7$']);%,num2str(hypers(2,1,1),'%0.1e'),'$'])

% name = 'lowRegGalaxie';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

ligne = 365;
ligne_mm = max(find(alpha_mm < alpha(ligne)));
load('/espace/orieux/results/fillGalaxie/coadd','coaddGalaxie')

% clf
% plot(beta, sky(ligne,:,2), 'r')
% hold on
% plot(beta, xchap(ligne,:,2))
% grid on

% xlabel('$\beta$')
% title(['$\gxb = 10^7$'])%,num2str(hypers(2,1,1),'%0.1e'),'$'])

% name = 'slicelowRegGalaxie';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%
% %% DSP

% sigma = 135;
% [ALPHA BETA] = ndgrid([-(Nalpha - 1)/2:(Nalpha - 1)/2],[-(Nbeta - 1)/2:(Nbeta - 1)/2]);
% window = exp(-1/2 * (ALPHA.^2 + BETA.^2)/sigma^2);

% noisefree = real(uifft2(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./Hrond360(1)));

% deltaF = 0.003;
% %% rf is reduced frequency
% [dspTrue rf] = estimCircularPSD(ufft2(window.*sky(:,:,2)),deltaF);
% dspProposed = estimCircularPSD(ufft2(window.*xchap(:,:,2)),deltaF);
% dspNoisefree = estimCircularPSD(ufft2(window.*noisefree),deltaF);
% dspFill = estimCircularPSD(ufft2(window.*fillmap(:,:,2)/Hrond360(1)),deltaF);

% f1 = 2;
% f2 = min(find(rf/2 >= 0.04));

% clf
% semilogy(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
% hold on
% semilogy(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
% semilogy(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
% semilogy(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')
% xlim([0 rf(f2)/alpha_step])
% grid on

% xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
% ylabel('$|S|^2$')
% legend('vrai', 'mcr', 'co-add', 'conv')
% legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')

% title(['$\gxb = 10^7$'])

% name = 'psdLowRegGalaxie';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

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
ylabel('$\alpha$')
xlabel('$\beta$')

name = 'bestGalaxie';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
plot(beta, sky(ligne,:,2), 'r')
hold on
plot(beta, xchap(ligne,:,2))
hold on
plot(beta_mm, coaddGalaxie(ligne_mm,:,2)/Hrond360(1),'k')
grid on
legend('vrai', 'mcr', 'co-add')
legend('Location',[0.6 0.6 0.32 0.25])

xlim([beta(1) beta(end)])

xlabel('$\beta$')

name = 'sliceBestGalaxie';
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

xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
ylabel('$|S|^2$')
legend('vrai', 'mcr', 'co-add', 'conv')
legend('Location',[0.2 0.2 0.32 0.25])
%legend('boxoff')

name = 'psdBestGalaxie';
laprint(1,[figplace,'/',name], 'width',8, 'asonscreen', 'on')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

% name = [placemount,expname,'/gpac_',num2str(Ntasks)];
% load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

% clf
% imagesc(alpha, beta, xchap(:,:,2))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')
% %title(['$\gxb = ',num2str(hypers(2,1,1),'%0.1e'),'$'])
% title(['$\gxb = 10^{15}$'])

% name = 'tomuchRegGalaxie';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

% clf
% plot(beta, sky(ligne,:,2), 'r')
% hold on
% % plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
% % hold on
% plot(beta, xchap(ligne,:,2))
% grid on

% xlabel('$\beta$')
% %title(['$\gxb = ',num2str(hypers(2,1,1),'%0.1e'),'$'])
% title(['$\gxb = 10^{15}$'])

% name = 'sliceTomuchRegGalaxie';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%
% %% DSP

% dspProposed = estimCircularPSD(ufft2(window.*xchap(:,:,2)),deltaF);

% clf
% semilogy(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
% hold on
% semilogy(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
% semilogy(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
% semilogy(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')
% xlim([0 rf(f2)/alpha_step])
% grid on

% xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
% ylabel('$|S|^2$')
% legend('vrai', 'mcr', 'co-add', 'conv')
% legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')
% title(['$\gxb = 10^{15}$'])

% name = 'psdTomuchRegGalaxie';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%% Comp madmap

% name = [placemount,expname,'/gpac_',num2str(theBest)];
% load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

% name = [placemount,'fillGalaxie/optim500000'];
% load(name, 'fillmap');

% clf
% plot(beta, sky(ligne,:,2), 'r')
% hold on
% plot(beta, fillmap(ligne,:,2)/Hrond360(1),'k')
% hold on
% plot(beta, xchap(ligne,:,2))
% grid on
% legend('Vrai', 'Co-add', 'MCR')
% legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')

% xlabel('$\beta$')

% name = 'sliceCompCoaddDeconvGalaxie';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%
% %%

% disp('Variance')

% Ntasks = 35;

% mean = zeros(Nalpha, Nbeta);
% cumulant2 = zeros(Nalpha, Nbeta);

% for itask = 1:Ntasks
%     for iworker = 1:8
%         name = [placemount,'varianceGalaxieDCT','/sample_',num2str(itask),'_',num2str(iworker)];
%         load(name, 'skySample')
        
%         mean = mean + conv2(skySample(:,:,2),fgaussian,'same');
%         cumulant2 = cumulant2 + (conv2(skySample(:,:,2),fgaussian,'same')).^2;
%     end
% end

% mean = mean/(Ntasks*8);
% cumulant2 = cumulant2/(Ntasks*8);

% variance = cumulant2 - mean.^2;
% standartDev = sqrt(variance);

% clf
% imagesc(alpha,beta,variance)
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')
% title(['$\gxb = 3,9\times 10^{12}$'])

% name = 'variance360Galaxie';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])


%%

% clf
% imagesc(alpha_mm, beta_mm, madmap(:,:,2)/Hrond360(1))
% axis image
% axis xy
% colorbar
% colormap(gray)

% ylabel('$\alpha$')
% xlabel('$\beta$')

% name = 'coadd';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

% clf
% plot(beta, fillmap(ligne,:,2)/Hrond360(1), 'r')
% hold on
% plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
% grid on

% xlabel('$\beta$')

% name = 'sliceCompCoaddFill';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

close all
