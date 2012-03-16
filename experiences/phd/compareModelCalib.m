% Comparison of the reponse model and the calibration data.

clear all

% format long
path(path,'/home/orieux/matlab/spire/trunk/libspire')
path(path,'/home/orieux/matlab/spire/trunk/simulator')
path(path,'/home/orieux/matlab/spire/trunk/utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

paramsInstrument
paramsObservation; simulatePointing;
paramsSky

radius_step = 0.005;

%% Precomputation
[alpha, beta, pointing250, pointing360, pointing520, index250, index360, ...
 index520, coefs250, coefs360, coefs520] = computeIndex(alpha, beta, ...
                                                  pointing250, pointing360, ...
                                                  pointing520, alpha_step, ...
                                                  beta_step);

precalculsInvariantRI

slice_model_250 = estimCircularPSD(fftshift(Hdirect250(:,:,1,1)),radius_step);

slice_model_360 = estimCircularPSD(fftshift(Hdirect360(:,:,1,1)),radius_step);

slice_model_520 = estimCircularPSD(fftshift(Hdirect520(:,:,1,1)),radius_step);

% Calibration
sigma_calib_psw = 18.1/sqrt(log(256));
sigma_calib_pmw = 25.2/sqrt(log(256));
sigma_calib_plw = 36.9/sqrt(log(256));

[ALPHAS BETAS] = ndgrid(SupAlpha, SupBeta);

psf_psw = 1/(2*pi*sigma_calib_psw^2)*exp(-(ALPHAS.^2 + BETAS.^2)./(2* ...
                                                  sigma_calib_psw^2));

psf_pmw = 1/(2*pi*sigma_calib_pmw^2)*exp(-(ALPHAS.^2 + BETAS.^2)./(2* ...
                                                  sigma_calib_pmw^2));

psf_plw = 1/(2*pi*sigma_calib_plw^2)*exp(-(ALPHAS.^2 + BETAS.^2)./(2* ...
                                                  sigma_calib_plw^2));

slice_calib_250 = estimCircularPSD(fftshift(psf_psw),radius_step);

slice_calib_360 = estimCircularPSD(fftshift(psf_pmw),radius_step);

slice_calib_520 = estimCircularPSD(fftshift(psf_plw),radius_step);

r = linspace(0,sqrt(2)*max(SupAlpha),size(slice_calib_250,1));

figure(1)
clf

subplot(4,3,1)
imagesc(Hdirect250(:,:,1,1))
axis square 
xlabel('\alpha')
ylabel('\beta')
title('PSW (model)')

subplot(4,3,2)
imagesc(Hdirect360(:,:,1,1))
axis square
xlabel('\alpha')
ylabel('\beta')
title('PMW (model)')

subplot(4,3,3)
imagesc(Hdirect520(:,:,1,1))
axis square
xlabel('\alpha')
ylabel('\beta')
title('PLW (model)')

subplot(4,3,4)
imagesc(psf_psw)
axis square
xlabel('\alpha')
ylabel('\beta')
title('PSW (calib)')

subplot(4,3,5)
imagesc(psf_pmw)
axis square
xlabel('\alpha')
ylabel('\beta')
title('PMW (calib)')

subplot(4,3,6)
imagesc(psf_plw)
axis square
xlabel('\alpha')
ylabel('\beta')
title('PLW (calib)')

% subplot(4,3,7)
% plot(r,slice_model_250)
% hold on
% plot(r,slice_calib_250./max(slice_calib_250).*max(slice_model_250),'--r')
% legend('Model','Calib')
% title('PSW Circular Mean')

% subplot(4,3,8)
% plot(r,slice_model_360)
% hold on
% plot(r,slice_calib_360./max(slice_calib_360).*max(slice_model_360),'--r')
% legend('Model','Calib')
% title('PSW Circular Mean')

% subplot(4,3,9)
% plot(r,slice_model_520)
% hold on
% plot(r,slice_calib_520./max(slice_calib_520).*max(slice_model_520),'--r')
% legend('Model','Calib')
% title('PSW Circular Mean')

slice_model_250 = Hdirect250(:,45,1,1);
slice_model_360 = Hdirect360(:,45,1,1);
slice_model_520 = Hdirect520(:,45,1,1);

slice_calib_250 = psf_psw(:,45,1,1);
slice_calib_360 = psf_pmw(:,45,1,1);
slice_calib_520 = psf_plw(:,45,1,1);

subplot(4,3,7)
plot(SupAlpha,slice_model_250)
hold on
plot(SupAlpha,slice_calib_250./max(slice_calib_250).*max(slice_model_250),'--r')
legend('Model','Calib')
xlabel('\beta (\alpha=0)')

subplot(4,3,8)
plot(SupAlpha,slice_model_360)
hold on
plot(SupAlpha,slice_calib_360./max(slice_calib_360).*max(slice_model_360),'--r')
legend('Model','Calib')
xlabel('\beta (\alpha=0)')

subplot(4,3,9)
plot(SupAlpha,slice_model_520)
hold on
plot(SupAlpha,slice_calib_520./max(slice_calib_520).*max(slice_model_520),'--r')
legend('Model','Calib')
xlabel('\beta (\alpha=0)')

%%


slice_model_250 = Hdirect250(45,:,1,1);
slice_model_360 = Hdirect360(45,:,1,1);
slice_model_520 = Hdirect520(45,:,1,1);

slice_calib_250 = psf_psw(45,:,1,1);
slice_calib_360 = psf_pmw(45,:,1,1);
slice_calib_520 = psf_plw(45,:,1,1);

subplot(4,3,10)
plot(SupAlpha,slice_model_250)
hold on
plot(SupAlpha,slice_calib_250./max(slice_calib_250).*max(slice_model_250),'--r')
legend('Model','Calib')
xlabel('\alpha (\beta=0)')

subplot(4,3,11)
plot(SupAlpha,slice_model_360)
hold on
plot(SupAlpha,slice_calib_360./max(slice_calib_360).*max(slice_model_360),'--r')
legend('Model','Calib')
xlabel('\alpha (\beta=0)')

subplot(4,3,12)
plot(SupAlpha,slice_model_520)
hold on
plot(SupAlpha,slice_calib_520./max(slice_calib_520).*max(slice_model_520),'--r')
legend('Model','Calib')
xlabel('\alpha (\beta=0)')
