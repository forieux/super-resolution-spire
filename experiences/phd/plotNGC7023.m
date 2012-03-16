placemount = '/espace/orieux/results/'
expname = 'ngc7023-2'

% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

figplace = '/home/orieux/spire/realData';
system(['mkdir -p ',figplace]);

set(0,'defaulttextinterpreter','tex')

set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxeslinewidth',1);
set(0,'defaultlinelinewidth',2);

SizeMarker = 2;
SizeLabel = 30;

load([placemount,'ngc7023-2/initialization'],'coaddCirrus2', 'alpha_mm', ...
     'beta_mm', 'alpha', 'beta', 'Hrond250', 'Hrond360', 'Hrond520', ...
     'alpha_step');

ligne = 840;
ligne_mm = min(find(alpha_mm > alpha(ligne)));

% for itask = 1:16
%     disp(['Run : ',num2str(itask)])
%     figure(1)

%     load([placemount,expname,'/gpac_PSW',num2str(itask)])
%     subplot(4,4,itask)
%     imagesc(log(xchap(:,:,1)+2)); axis equal

%     figure(2)
%     load([placemount,expname,'/gpac_PMW',num2str(itask)])
%     subplot(4,4,itask)
%     imagesc(log(xchap(:,:,2)+2)); axis equal

%     figure(3)
%     load([placemount,expname,'/gpac_PLW',num2str(itask)])
%     subplot(4,4,itask)
%     imagesc(log(xchap(:,:,3)+2)); axis equal
    
% %     drawnow
% %     pause
    
% %     figure(2)
% %     clf
% %     plot(beta, xchap(ligne,:,1),'r');
% %     hold on
% %     plot(beta_mm, coaddCirrus2(ligne_mm,:,1)/Hrond250(1))    
% end

load([placemount,expname,'/gpac_PLW12'])
a = xchap(:,:,3);
save PLW a -ASCII

load([placemount,expname,'/gpac_PMW11'])
a = xchap(:,:,2);
save PMW a -ASCII

load([placemount,expname,'/gpac_PSW10'])
a = xchap(:,:,1);
save PSW a -ASCII


break

load /espace/orieux/results/ngc7023-2/gpac_PLW12.mat
PLW = xchap(:,:,3);
load /espace/orieux/results/ngc7023-2/gpac_PMW11.mat
PMW = xchap(:,:,2);
load /espace/orieux/results/ngc7023-2/gpac_PSW9.mat
PSW = xchap(:,:,1);

% PLW=PLW-min2(PLW)+min2(PSW);
% PMW=PMW-min2(PMW)+min2(PSW);
% PSW=PSW./max2(PLW)*max2(PSW);

PLW=PLW.*Hrond520(1);
PMW=PMW.*Hrond360(1);
PSW=PSW.*Hrond250(1);

ligne=935;
ligne_mm = min(find(alpha_mm > alpha(ligne)));

figure(3)
clf
plot(beta,PSW(ligne,:),'b')
hold on
plot(beta,PMW(ligne,:),'g')
plot(beta,PLW(ligne,:),'r')
plot(beta_mm,coaddCirrus2(ligne_mm,:,1),'b--')
hold on
plot(beta_mm,coaddCirrus2(ligne_mm,:,2),'g--')
plot(beta_mm,coaddCirrus2(ligne_mm,:,3),'r--')
xlim([beta(400) beta(1400)])
legend('PSW', 'PMW', 'PLW', 'PSW coadd', 'PMW coadd', 'PLW coadd')

save PSW PSW -ASCII
save PMW PMW -ASCII
save PLW PLW -ASCII
save alpha alpha -ASCII
save beta beta -ASCII

break

% %%
% %% DSP

boundAlpha = [450 1350] ;
boundBeta = [400 1300];

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

deltaF = 0.004;


%% rf is reduced frequency

load /espace/orieux/results/ngc7023/gpac_7
proposed = xchap(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2),1)*Hrond250(1);

load /espace/orieux/results/ngc7023/optim10000
filled = fillmap(boundAlpha(1):boundAlpha(2), boundBeta(1):boundBeta(2),1);
% filled = filled./mean2(filled)*mean2(proposed);

[dspProposed rf] = estimCircularPSD(ufft2(window.*proposed),deltaF);
dspFill = estimCircularPSD(ufft2(window.*filled),deltaF);

f1 = min(find(rf/2 >= 0.001));
f2 = min(find(rf/2 >= 0.15));

clf
% RF = F/FE = F*Te

loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
hold on
loglog(rf(f1:f2)/alpha_step, dspFill(f1:f2),'k')
xlim([0.001 0.11])
ylim([1e-5 1e4])
% % xlim([0.005 0.11])
% % ylim([1e-17 1e-7])
grid on

xlabel('1/arcsec')
ylabel('|S|^2', 'interpreter','tex')
legend('proposed', 'co-add')

% name = 'dsp'
% print(1,'-depsc',[figplace,'/',name])

%%

% proposed = xchap(:,:,1);
% save proposed proposed -ASCII -DOUBLE
% coadd = coaddCirrus2(:,:,1);
% save coadd coadd -ASCII -DOUBLE

%%

dynamique = [min2(xchap(:,:,1)) max2(xchap(:,:,1))];
figure(1)
clf
imagesc(alpha, beta, xchap(:,:,1), dynamique)
xlim([3.928 3.9435]*1e5)
ylim([2.4468 2.461]*1e5)
axis off
axis equal

name = 'imProposed';
print(1,'-depsc',[figplace,'/',name])

figure(1)
clf
imagesc(alpha_mm, beta_mm, coaddCirrus2(:,:,1)/Hrond250(1), dynamique)
xlim([3.928 3.9435]*1e5)
ylim([2.4468 2.461]*1e5)
axis off
axis equal

name = 'imCoadd';
print(1,'-depsc',[figplace,'/',name])

%%

alphaBox = [2.451 2.4555]*1e5;
betaBox = [3.934 3.94]*1e5;

figure(1)
clf
imagesc(alpha, beta, xchap(:,:,1), dynamique)
ylim([2.451 2.4555]*1e5)
xlim([3.934 3.94]*1e5)
axis off
axis equal

name = 'zoomProposed'
print(1,'-depsc',[figplace,'/',name])

figure(1)
clf
imagesc(alpha_mm, beta_mm, coaddCirrus2(:,:,1)/Hrond250(1), dynamique)
ylim([2.451 2.4555]*1e5)
xlim([3.934 3.94]*1e5)
axis off
axis equal

name = 'zoomCoadd'
print(1,'-depsc',[figplace,'/',name])

%%

ligne = 840;
ligne_mm = min(find(alpha_mm > alpha(ligne)));

figure(1)
clf
plot(beta, xchap(ligne,:,1),'r')
hold on
plot(beta_mm, coaddCirrus2(ligne_mm,:,1)/Hrond250(1))
xlim([2.445 2.46]*1e5)
grid on

legend('proposed', 'co-add')

name = 'slice'
print(1,'-depsc',[figplace,'/',name])
