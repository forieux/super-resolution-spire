% format long

path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

figplace = '/home/orieux/tex/manuscrit/imagesRes/prior/';
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

% load('/espace/orieux/results/variousDirty/dirtymaps')
% load('/espace/orieux/results/simpleFillCirrus/optim50000','fillmap')

placemount = '/espace/orieux/results/'
expname = 'variousRegPriorDCT';
Ntasks = 144;

boundAlpha = [160 420];
boundBeta = [160 420];

boundBetal = [160 420];

load([placemount,expname,'/initialization'])

%%

load /espace/orieux/results/bestRegPrior/gpac

l1PriorEAP = sum2(abs(sky(:,:,2) - xchap(:,:,2)))/sum2(abs(sky(:,:,2))) 
l2PriorEAP = sum2((sky(:,:,2) - xchap(:,:,2)).^2)/sum2((sky(:,:,2).^2)) 

%%

true = conv2(sky(:,:,2),fgaussian,'same');

clf
imagesc(alpha, beta, true)
axis image
axis xy
colorbar
colormap(gray)
ylabel('$\beta$')
xlabel('$\alpha$')

name = 'imagePrior';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%


% clf
% imagesc(alpha, beta, mapmean(:,:,1)/Hrond250(1))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')

% name = 'dirtyMean250';
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

% name = 'dirtyMean360';
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

% name = 'dirtyMean520';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

% clf
% imagesc(alpha, beta, fillmap(:,:,1)/Hrond250(1))
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')

% name = 'fillMap250';
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

% name = 'fillMap360';
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

% name = 'fillMap520';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%% L2 et L1
for itask = 1:Ntasks
    
    disp(num2str(itask));
    
    name = [placemount,expname,'/gpac_',num2str(itask)];
    load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')
    
    dist1(itask) = sum2(abs(conv2(sky(:,:,2) - xchap(:,:,2),fgaussian, ...
                                  'same')))/sum2(abs(conv2(sky(:,:,2), ...
                                                      fgaussian,'same'))); 
    dist2(itask) = sum2((conv2(sky(:,:,2) - xchap(:,:,2),fgaussian, ...
                               'same')).^2)/sum2((conv2(sky(:,:,2), ...
                                                      fgaussian,'same').^2)); 
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


% clf
% semilogx(sorted1(:,1),sorted1(:,2))
% hold on
% semilogx(sorted2(:,1),sorted2(:,2),'r')
% grid on
% legend('$\ell_1$','$\ell_2$')
% legend('boxoff')
% plot(sorted1(theBest1,1), sorted1(theBest1,2),'.');%, 'MarkerSize', SizeMarker)
% plot(sorted2(theBest2,1), sorted2(theBest2,2),'.r');%, 'MarkerSize', SizeMarker)

% xlim([min(sorted1(:,1)) max(sorted1(:,1))])

% xlabel('$\gxb$')

% name = 'bestHyperPrior';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

% name = [placemount,expname,'/gpac_',num2str(1)];
% load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')
% sol = conv2(xchap(:,:,2),fgaussian,'same');
% true = conv2(sky(:,:,2),fgaussian,'same');

% clf
% imagesc(alpha, beta, sol)
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')
% title(['$\gxb = 10^7$']);%,num2str(hypers(2,1,1),'%0.1e'),'$'])

% name = 'lowRegPrior';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

% ligne = 200;
% load('/espace/orieux/results/fillPrior/coadd','coaddPrior','alpha_mm', 'beta_mm')
% ligne_mm = max(find(alpha_mm < alpha(ligne)));

% clf
% plot(beta, true(ligne,:), 'r')
% hold on
% plot(beta, sol(ligne,:))

% grid on

% xlabel('$\beta$')
% title(['$\gxb = 10^7$'])%,num2str(hypers(2,1,1),'%0.1e'),'$'])

% name = 'slicelowRegPrior';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% DSP

% sigma = 135;
% [ALPHA BETA] = ndgrid([-(Nalpha - 1)/2:(Nalpha - 1)/2],[-(Nbeta - 1)/2:(Nbeta - 1)/2]);
% window = exp(-1/2 * (ALPHA.^2 + BETA.^2)/sigma^2);

% noisefree = real(uifft2(ufft2(sky(:,:,2).*Hrond360(:,:,1)./Hrond360(1))));

% deltaF = 0.003;
% %% rf is reduced frequency
% [dspTrue rf] = estimCircularPSD(ufft2(window.*true),deltaF);
% dspProposed = estimCircularPSD(ufft2(window.*sol),deltaF);
% dspNoisefree = estimCircularPSD(ufft2(window.*noisefree),deltaF);

% f1 = 2;
% f2 = min(find(rf/2 >= 0.04));

% clf
% semilogy(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
% hold on
% semilogy(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
% semilogy(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')
% xlim([0 rf(f2)/alpha_step])
% grid on

% xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
% ylabel('$|S|^2$')
% legend('vrai', 'mcr', 'co-add', 'conv')
% legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')

% title(['$\gxb = 10^7$'])

% name = 'psdLowRegPrior';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])


%%
%%

lesRegs = logspace(7,15,144);
theBest = find(lesRegs == sorted1(find(sorted1(:,2) == min(sorted1(:,2))),1));

name = [placemount,expname,'/gpac_',num2str(theBest)];
load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')
sol = conv2(xchap(:,:,2),fgaussian,'same');

clf
imagesc(alpha, beta, sol)
axis image
axis xy; axis off
colorbar
colormap(gray)

title(['$\gxb = 7\times 10^{11}$'])

name = 'bestPrior';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

% clf
% plot(beta, true(ligne,:), 'r')
% hold on
% plot(beta, sol(ligne,:))
% grid on

% xlabel('$\beta$')
% title(['$\gxb = 7\times 10^{11}$'])

% name = 'sliceBestPrior';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% DSP

% dspProposed = estimCircularPSD(ufft2(window.*sol),deltaF);

% clf
% semilogy(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
% hold on
% semilogy(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
% %semilogy(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
% semilogy(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')
% xlim([0 rf(f2)/alpha_step])
% grid on

% xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
% ylabel('$|S|^2$')
% legend('vrai', 'mcr', 'co-add', 'conv')
% legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')
% title(['$\gxb = 1,4\times 10^{12}$'])

% name = 'psdBestPrior';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

% name = [placemount,expname,'/gpac_',num2str(Ntasks)];
% load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')
% sol = conv2(xchap(:,:,2),fgaussian,'same');

% clf
% imagesc(alpha, beta, sol)
% axis image
% axis xy
% colorbar
% colormap(gray)
% ylabel('$\alpha$')
% xlabel('$\beta$')
% %title(['$\gxb = ',num2str(hypers(2,1,1),'%0.1e'),'$'])
% title(['$\gxb = 10^{15}$'])

% name = 'tomuchRegPrior';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

% clf
% plot(beta, true(ligne,:), 'r')
% hold on
% plot(beta, sol(ligne,:))
% grid on

% xlabel('$\beta$')
% %title(['$\gxb = ',num2str(hypers(2,1,1),'%0.1e'),'$'])
% title(['$\gxb = 10^{15}$'])

% name = 'sliceTomuchRegPrior';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% DSP

% dspProposed = estimCircularPSD(ufft2(window.*sol),deltaF);

% clf
% semilogy(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
% hold on
% semilogy(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
% %semilogy(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
% semilogy(rf(f1:f2)/alpha_step, dspNoisefree(f1:f2),'--k')
% xlim([0 rf(f2)/alpha_step])
% grid on

% xlabel('$f~\big[\textrm{arcsec}^{-1}\big]$')
% ylabel('$|S|^2$')
% legend('vrai', 'mcr', 'co-add', 'conv')
% legend('Location',[0.2 0.2 0.32 0.25])
% legend('boxoff')
% title(['$\gxb = 10^{15}$'])

% name = 'psdTomuchRegPrior';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

% disp('Variance')

% Ntasks = 35;

% mean = zeros(Nalpha, Nbeta);
% cumulant2 = zeros(Nalpha, Nbeta);

% for itask = 1:Ntasks
%     for iworker = 1:8
%         name = [placemount,'variancePriorDCT','/sample_',num2str(itask),'_',num2str(iworker)];
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
% title(['$\gxb = 7\times 10^{11}$'])

% name = 'variance360Prior';
% laprint(1,[figplace,'/',name], 'width',10)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% Correlation coef

%%
%%
%% Comp madmap

load([placemount,expname,'/gpac_88'],'xchap')
sol = conv2(xchap(:,:,2),fgaussian,'same');

ligne = 200;
load('/espace/orieux/results/fillPrior/coadd','coaddPrior','alpha_mm', 'beta_mm')
ligne_mm = max(find(alpha_mm < alpha(ligne)));

clf
plot(beta(10:Nbeta), true(ligne,10:Nbeta), 'r')
hold on
plot(beta_mm, coaddPrior(ligne_mm,:,2)/Hrond360(1),'k')
plot(beta(10:Nbeta), sol(ligne,10:Nbeta))
grid on

xlim([beta(10) beta(end)])
ylim([2 7]*1e-5)

xlabel('$\beta$')
title(['$\gxb = 7^{11}$'])

name = 'sliceCompCoaddDeconvPrior';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])


% name = [placemount,expname,'/gpac_',num2str(theBest)];
% load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')

% clf
% plot(beta, sky(ligne,:,2), 'r')
% hold on
% plot(beta_mm, madmap(ligne_mm,:,2)/Hrond360(1),'k')
% hold on
% plot(beta, xchap(ligne,:,2))
% grid on

% xlabel('$\beta$')

% name = 'sliceCompCoaddDeconv';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% %%

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

% %%

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
