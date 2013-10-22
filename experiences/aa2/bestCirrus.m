clear all

%% Base
placemount = '/mnt/space/results/';
expname = 'aa2-bestCirrus';
savedir = [placemount, expname, '/']
system(['mkdir -p ',savedir]);
addpath('../../')
addpath('../../utils')
addpath('../../simulator')
addpath('../../libspire')

randn('state',0)
rand('state',0)

%%
paramsInstrument
paramsObservation;
simulatePointing;
paramsSky

%% Precomputation
bound = 4*sigma_coef*max(band_520) + 60; % arcsec
[alpha beta] = computeAxis(pointing250, sky_alpha_period, ...
                           sky_beta_period, ...
                           bound, N_scan_total);

%% Precomputation of index and redondancy
[alphaCoadd, betaCoadd, pointing250coadd, index250coadd, ...
 coefs250coadd, Nalphacoadd, Nbetacoadd] = computeIndex(alpha, beta, ...
                                                        pointing250, ...
                                                        6, 6, ...
                                                        Nbolo250, ...
                                                        N_scan_total);

[alpha, beta, pointing250, index250, coefs250, Nalpha, Nbeta] = ...
    computeIndex(alpha, beta, pointing250, alpha_step, beta_step, ...
                 Nbolo250, N_scan_total);

%% The coaddition without offsets and estimation of offsets
offsets = zeros(1,Nbolo250);

%% Inversion
init = zeros(Nalpha, Nbeta);
objectMean = zeros(size(init));

[Hrond250 Hdirect250] = computeRI(10, 1, band_250, ...
                                  central_wavelength_250, ...
                                  sigma_coef, sigma_alpha, ...
                                  sigma_beta, alpha_step, beta_step, ...
                                  Nalpha, Nbeta, Nspeed, Norder, ...
                                  unique_speed, opt_efficiency, ...
                                  time_constante, gain);
%% Gain
G = Hrond250(1);

[diffAlpha diffBeta circDalpha circDbeta] = computeReg(sigma_alpha, ...
                                                       sigma_beta, ...
                                                       sky_alpha_period, ...
                                                       sky_beta_period, ...
                                                       5, Nalpha, ...
                                                       Nbeta);

regOp = circDalpha + circDbeta;
regOp(1) = 0;

%% Simulate data
std250 = (1e-3); gammaB250 = 1/std250^2;
std360 = std250; gammaB360 = 1/std360^2;
std520 = std360; gammaB520 = 1/std520^2;

std = [std250 std360 std520];
gammaB = [gammaB250, gammaB360, gammaB520];

[data250 output sky] = ...
    simulateData(['../../data&true/simul_ism_SPS_comp'], Nalpha, ...
                 Nbeta, Norder, std250, 10^(-4), Hrond250, ...
                 index250, Nbolo250, Nspeed, N_scan_total, unique_speed, the_speeds);

coaddPSW = dirtymap(data250, index250coadd, zeros(Nbolo250, 1), Nalphacoadd, Nbetacoadd, N_scan_total, Nbolo250);
coadd_interp = imresize(coaddPSW, size(sky), 'nearest');
mesure_sky = sky;
mesure_sky(find(coadd_interp == 0)) = 0;
disp(num2str(sum(abs(coadd_interp(:)/Hrond250(1) - mesure_sky(:)).^2) / sum(abs(mesure_sky(:)).^2)))

%% Options
cgoptions.thresold = 1e-8;
cgoptions.maxIter = 40;
cgoptions.numfig = 1000;
%%
init = 6.5598e-05 * ones(Nalpha, Nbeta);

objectMean = zeros(size(init));
gamma_x = logspace(11, 14, 100);
error = zeros(size(gamma_x));
best_error = Inf;
for iIter = 1:length(gamma_x)
  hypers = [gammaB250; gamma_x(iIter)];
  [skyMapF ooptim] = inversionF(ufft2(init), cgoptions, data250, hypers, ...
                                  Hrond250, index250, coefs250, offsets, ...
                                  regOp, objectMean, Nalpha, Nbeta, ...
                                  Norder, N_scan_total, Nbolo250, Nspeed, ...
                                  unique_speed, the_speeds);

  skyMap =  real(uifft2(skyMapF));
  imagesc(skyMap)
  mesure_skyMap = skyMap;
  mesure_skyMap(find(coadd_interp == 0)) = 0;

  error(iIter) = sum(abs(mesure_skyMap(:) - mesure_sky(:)).^2) / sum(abs(mesure_sky(:)).^2);
  if error(iIter) < best_error
    best_error = error(iIter);
    optimSky = skyMap;
  end
end

%%
[ALPHAS BETAS] = ndgrid([-5:5], [-5:5]);
sigma_alpha = 0.6*alpha_step;
sigma_beta = 0.6 * beta_step;
fgaussian = 1/(2*pi*sigma_alpha * sigma_beta)* exp(- ALPHAS.^2 / (2 * sigma_alpha^2) ...
                                           - BETAS.^2 / (2* sigma_beta^2));

%% Unsupervised estimate
criterion = 1e-4;
burnin = 1000; %% If you don't know this value it is more than 100
maxIter = 2500; %% If you don't know it is more 200

hypersInit = zeros(2,1);
hypersInit(1,1) = 1e5;
hypersInit(2,1) = 1e9;

[skyEap gnChain gxChain skyStd] = usmse(init, hypersInit, data250, Hrond250, index250, coefs250, offsets, regOp, Nalpha, Nbeta, Norder, N_scan_total, Nbolo250, Nspeed, unique_speed, the_speeds, criterion, burnin, maxIter, cgoptions, fgaussian, 'pla', [placemount,expname]);

save([savedir, 'workspace.mat'])

break

figplace = '/home/orieux/tex/papers/AA2/tex/figs2';
system(['mkdir -p ',figplace]);

optim_disp = conv2(optimSky, fgaussian, 'same');
eap_disp = conv2(skyEap, fgaussian, 'same');
sky_disp = conv2(sky, fgaussian, 'same');

set(0,'defaultaxesfontsize',25);
set(0,'defaulttextfontsize',25);
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);

SizeMarker = 2;
SizeLabel = 30;

numfig = 1;

%%
sc = [min(skyEap(:)) max(skyEap(:))];

figure(2)
subplot(121)
imagesc(skyEap, sc)
axis image
axis off
colormap gray
subplot(122)
imagesc(optimSky, sc)
axis image
colormap gray
axis off

figure(numfig)
clf()
sc = [min(eap_disp(:)) max(eap_disp(:))];
imagesc(eap_disp, sc)
axis image
axis off
colormap gray
print('-depsc',[figplace,'/eapCirrus']);

clf
imagesc(optim_disp, sc)
axis image
colormap gray
axis off
print('-depsc',[figplace,'/bestCirrus']);

clf
imagesc(coadd_interp / Hrond250(1), sc)
axis image
colormap gray
axis off
print('-depsc',[figplace,'/coaddCirrus']);

xmin=100;
xmax=400;
line = 200;

clf()
plot(eap_disp(line,xmin:xmax))
hold on
plot(optim_disp(line,xmin:xmax), 'r')
plot(sky_disp(line,xmin:xmax), 'k')
grid on
legend('Proposed', 'Best', 'True')
xlim([0 (xmax - xmin)])
ylim([4 8] * 1e-5)
print('-depsc',[figplace,'/sliceCirrusMoy']);

clf()
plot(sky_disp(line,xmin:xmax))
hold on
plot(coadd_interp(line,xmin:xmax) / Hrond250(1), 'r')
legend('Proposed', 'Naive')
grid on
xlim([0 (xmax - xmin)])
ylim([4 8] * 1e-5)
print('-depsc',[figplace,'/sliceCirrusMoyNaive']);

clf()
sky_disp(line:line+2,xmin:xmax) = max(sky(:));
imagesc(sky_disp, sc)
axis image
colormap gray
axis off
print('-depsc',[figplace,'/cirrus']);


%%
% Window Gaussienne
sigma = 20;

boundAlpha = [160 420];
boundBeta = [160 420];
NalphaW = boundAlpha(2) - boundAlpha(1);
NbetaW = boundBeta(2) - boundBeta(1);

[ALPHA BETA] = ndgrid([-(NalphaW)/2:(NalphaW)/2],[-(NbetaW)/2:(NbetaW)/2]);
window = exp(-1/2 * (ALPHA.^2 + BETA.^2)/sigma^2);
%window = ones(size(window));

deltaF = 0.004;
[dspTrue rf] = estimCircularPSD(ufft2(window.*sky_disp(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);
dspProposed = estimCircularPSD(ufft2(window.*eap_disp(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);
dspOptim = estimCircularPSD(ufft2(window.*optim_disp(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);
dspCoadd = estimCircularPSD(ufft2(window.*coadd_interp(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))/Hrond250(1)),deltaF);

f1 = min(find(rf/2 >= 0.001));
f2 = min(find(rf/2 >= 0.15));

clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2))
loglog(rf(f1:f2)/alpha_step, dspOptim(f1:f2),'k')
loglog(rf(f1:f2)/alpha_step, dspCoadd(f1:f2),'--k')

xlim([0.005 0.11])
ylim([1e-17 1e-7])

grid on

xlabel('f [1/arcsec]')
legend('True', 'Proposed', 'Best', 'Naive')
legend('Location',[0.20 0.20 0.2 0.25])

print('-depsc',[figplace,'/psdEapCirrus'])

%%
clf()
imagesc(skyStd)
axis image
colormap gray
axis off
print('-depsc',[figplace,'/varianceCirrus']);

clf()
sc = [min(min(skyStd(200:400))) max(max(skyStd(200:400)))];
imagesc(abs(optim_disp(200:400,200:400) - sky_disp(200:400,200:400)))
imagesc(abs(optimSky(200:400,200:400) - sky(200:400,200:400)))
axis image
colormap gray
colorbar
axis off
print('-depsc',[figplace,'/optim_residual']);

clf()
imagesc(abs(skyEap(200:400,200:400) - sky(200:400,200:400)))
axis image
colormap gray
colorbar
axis off
print('-depsc',[figplace,'/eap_residual']);

clf()
imagesc(abs(coadd_interp(200:400,200:400)/Hrond250(1) - sky(200:400,200:400)))
axis image
colormap gray
colorbar
axis off
print('-depsc',[figplace,'/coadd_residual']);

break

% clf()
% plot(skyStd(line,:))
% print('-depsc',[figplace,'/sliceStd']);

clf()
plot(sky(line,xmin:xmax))
hold on
plot(eap_disp(line,xmin:xmax) - skyStd(line,xmin:xmax), 'k--')
plot(eap_disp(line,xmin:xmax) + skyStd(line,xmin:xmax), 'k--')
grid on
xlim([0 (xmax - xmin)])
ylim([4 8] * 1e-5)
print('-depsc',[figplace,'/sliceStandartDevGibbsCirrus']);

%%
mesure_sky = sky;
mesure_sky(coadd_interp == 0) = 0;
mesure_optim = optimSky;
mesure_optim(coadd_interp == 0) = 0;
disp(['Optim ', num2str(sum(abs(mesure_optim(:) - mesure_sky(:)).^2) / sum(abs(mesure_sky(:)).^2))])

mesure_eap = skyEap;
mesure_eap(coadd_interp == 0) = 0;
disp(['Eap ', num2str(sum(abs(mesure_eap(:) - mesure_sky(:)).^2) / sum(abs(mesure_sky(:)).^2))])

disp(['Naive ', num2str(sum(abs(coadd_interp(:)/Hrond250(1) - mesure_sky(:)).^2) / sum(abs(mesure_sky(:)).^2))])

%%
% iter = [1:length(gnChain(:))];
% currentMean = cumsum(gnChain(:));
%
% clf
% plot(gnChain(:))
% hold on
% xlim([iter(1)-100 iter(end)+100])
%
% grid on
%
% print('-depsc',[figplace,'/chaineGammaB-jfg'])
%
% clf
% hist(gnChain(400:end),50)
% xlim([0.97 1.03]*1e6)
%
% print('-depsc',[figplace,'/histGammaB-jfg'])
%
% meangn = mean(gnChain(400:end))
% stdgn = sqrt(mean(gnChain(400:end).^2) - meangn^2)

clf
subplot(1,2,2, 'Position', [0.65, 0.1, 0.25, .82])
hist(gnChain(:), linspace(9e5, 11e5, 40))
set(gca,'CameraUpVector', [1,0,0]);
xlim([-0 12e5])
set(gca, 'ytick', []);
grid on
subplot(1,2,1, 'Position', [0.1, 0.1, 0.6, .82])
plot(gnChain(:), '.')
hold on
%plot(cumsum(gnChain,:)) ./ cumsum(ones(size(offsetsChain(54,:)))), 'r')
ylim([-0 12e5])
xlim([-100 length(gnChain) + 100])
grid on
print('-depsc',[figplace,'/chaineGammaB-jfg'])

%%
iter = [1:length(squeeze(gxChain(:)))]';
% currentMean = cumsum(gxChain(500:end))./(iter(500:end) - 500);

% clf
% plot(squeeze(gxChain(:)))
% hold on
% xlim([iter(1) - 100 iter(end) + 100])
% ylim([min(gxChain(:))-1e11 max(gxChain(:))+1e11])
% grid on
%
% print('-depsc',[figplace,'/chaineGammaX-jfg'])
%
% clf
% hist(gxChain(600:end),50)
% xlim([8 10]*1e11)
%
% print('-depsc',[figplace,'/histGammaX-jfg'])
%
% disp(['Gx mean ', num2str(mean(gxChain(burnin:end)))])
% disp(['Gx std ', num2str(sqrt(mean(gxChain(burnin:end).^2) - mean(gxChain(burnin:end))^2))])
% %disp(['Best Gx ', num2str(gamma_x(error == best_error))])


clf
subplot(1,2,2, 'Position', [0.65, 0.1, 0.25, .82])
hist(gxChain(burnin:end), linspace(8e11, 10e11, 40))
set(gca,'CameraUpVector', [1,0,0]);
xlim([-1e11 11e11])
set(gca, 'ytick', []);
grid on
subplot(1,2,1, 'Position', [0.1, 0.1, 0.6, .82])
plot(gxChain(:), '.')
%hold on
%plot(cumsum(gnChain,:)) ./ cumsum(ones(size(offsetsChain(54,:)))), 'r')
ylim([-1e11 11e11])
xlim([-100 length(gxChain) + 100])
grid on
print('-depsc',[figplace,'/chaineGammaX-jfg'])

%%
name = [placemount,expname,'/sample_',num2str(1)];

[dspTrue rf] = estimCircularPSD(ufft2(window.*sky_disp(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);
dspProposed = estimCircularPSD(ufft2(window.*eap_disp(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2))),deltaF);

cumulant1r = zeros(size(dspProposed));
cumulant2r = zeros(size(cumulant1r));

shist = [];
fhist = [];
idx_min = find(log(rf(2:end)/alpha_step) >= log(0.005), 1, 'first');
idx_max = find(log(rf(2:end)/alpha_step) <= log(0.11), 1, 'last');
lf = linspace(log(rf(idx_min)/alpha_step), log(rf(idx_max)/alpha_step), 200);
lstrue = interp1(log(rf(idx_min:idx_max)/alpha_step), log(dspTrue(idx_min:idx_max)), lf, 'linear');
N = 1500;
for i = 1000:1000+N
    disp(num2str(i))
    load([placemount,expname,'/samples/sample_',num2str(i)], 'skySample')

    im = conv2(real(uifft2(skySample)), fgaussian, 'same');
    s = estimCircularPSD(ufft2(window .* im(boundAlpha(1):boundAlpha(2), boundBeta(1):boundBeta(2))), deltaF);
    ls = interp1(log(rf(idx_min:idx_max)/alpha_step), log(s(idx_min:idx_max)), lf, 'linear');

    shist = cat(1, shist, ls(:));
    fhist = cat(1, fhist, lf(:));
end
shist = cat(1, shist, lstrue(:));
fhist = cat(1, fhist, lf(:));

%%
[N, C] = hist3(cat(2, fhist, shist), [100, 100]);
ltrue = interp1(lf, lstrue, C{1});
lcoadd = interp1(lf, interp1(log(rf(idx_min:idx_max)/alpha_step), log(dspCoadd(idx_min:idx_max)), lf, 'linear'), C{1});
lproposed = interp1(lf, interp1(log(rf(idx_min:idx_max)/alpha_step), log(dspProposed(idx_min:idx_max)), lf, 'linear'), C{1});
% lprop = interp1(lf, interp1(log(rf(idx_min:idx_max)/alpha_step), log(dspProposed(idx_min:idx_max)), lf, 'linear'));
% lstd = interp1(lf, interp1(log(rf(idx_min:idx_max)/alpha_step), log(standartDeviationR(idx_min:idx_max)), lf, 'linear'));

clf()
imagesc(C{1}, C{2}, log(N'))
hold on
plot(C{1}, ltrue, 'r')
hold on
plot(C{1}, lcoadd, 'k--')
grid on
% plot(C{1}, lproposed)
legend('True', 'Coadd')

% caxis([0, 0.5 * max(N(:))])
colormap(1-gray)
shading interp
axis xy;
xlim([min(C{1}) max(C{1})])
ylim([min(C{2}) max(C{2})])

xlim(log([0.005 0.11]))
ylim(log([1e-17 1e-7]))

set(gca,'XTickLabel', num2cell(num2str(reshape(exp(get(gca, 'Xtick')), [], 1),'%10.1e\n' ), 2))
set(gca,'YTickLabel', num2cell(num2str(reshape(exp(get(gca, 'Ytick')), [], 1),'%10.1e\n' ), 2))
xlabel('f [1/arcsec]')

print('-depsc',[figplace,'/psdCirrusDensity'])
%%

xtlabel = {}
for xt = reshape(exp(get(gca, 'Xtick')), [], 1)
    print xt
end

set(gca,'XTickLabel',['1 ';'100'])
%%
clf
loglog(rf(f1:f2)/alpha_step, dspTrue(f1:f2),'r')
hold on
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2) + 1*standartDeviationR(f1:f2),'k--')
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2) + 3*standartDeviationR(f1:f2),'k.')
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2) - 1*standartDeviationR(f1:f2),'k--')
loglog(rf(f1:f2)/alpha_step, dspProposed(f1:f2) - 3*standartDeviationR(f1:f2),'k.')
xlim([0.005 0.11])
ylim([1e-17 1e-7])
legend('True', '\pm \sigma', '\pm 3 \sigma')
grid on
xlabel('f [1/arcsec]')

print('-depsc',[figplace,'/psdCirrusStd'])

%%
hits = zeros(size(coefs250{1}));
for i = 1:length(coefs250)
    hits = hits + coefs250{i};
end
hits(1:2,:) = max(hits(:));
hits(end-1:end,:) = max(hits(:));
hits(:,1:2) = max(hits(:));
hits(:,end-1:end) = max(hits(:));

clf()
imagesc(hits)
axis image
colormap(1-gray)
colorbar
axis off
print('-depsc',[figplace,'/redondance250'])

%%


%%
