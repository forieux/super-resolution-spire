% Time-stamp: <2010-11-15 11:32:43 orieux>

alpha_step = 2;
beta_step = 2;

sky_alpha_period = alpha_step;
sky_beta_period  = sky_alpha_period;

%%% 0.6 is the best to approx the equivalent sinc function...
sigma_alpha = 0.6*sky_alpha_period;
sigma_beta  = 0.6*sky_beta_period;

%%% Alpha and beta coordinate
% extra space for the convolution + size of the sensor 4'x8'
bord = 4*sigma_coef*max(band_520) + 4*60; % in arcsec

p = pointing_band{1};
minAlpha = min(min(p(1,:,:)));
minBeta = min(min(p(2,:,:)));
maxAlpha = max(max(p(1,:,:)));
maxBeta = max(max(p(2,:,:)));
for iscan = 2:N_scan_total
    p = pointing_band{iscan};
    minAlpha = min(minAlpha, min(min(p(1,:,:))));
    minBeta  = min(minBeta,  min(min(p(2,:,:))));
    maxAlpha = max(maxAlpha, max(max(p(1,:,:))));
    maxBeta  = max(maxBeta,  max(max(p(2,:,:))));
end

clear p
alpha = [minAlpha - bord : sky_alpha_period : maxAlpha + bord];
beta = [minBeta - bord : sky_beta_period : maxBeta + bord];

%The RI will be computed on all the support of the sky, to avoid the phase
%problem in Fourrier space. Secondly to be certain that the RI will be
%computed on zero coordinate, we need odd number of element.
if mod(length(alpha),2) == 0
    alpha = [alpha max(alpha)+sky_alpha_period];
end
if mod(length(beta),2) == 0
    beta = [beta max(beta)+sky_beta_period];
end

Nalpha = length(alpha);
Nbeta = length(beta);

delta_alpha = alpha(2) - alpha(1);
delta_beta = beta(2) - beta(1);

%% The number of order (so 1 if you want only order 0)
Norder = 1;

%% Gaussian decomposition function
alphaS = [-5*sigma_alpha:sky_alpha_period:5*sigma_alpha];
[ALPHAS BETAS] = ndgrid(alphaS,alphaS);

%% Fonction de d�composition 
fgaussian = 1/(2*pi*sigma_alpha*sigma_beta)* exp(-ALPHAS.^2/(2*sigma_alpha^2) ...
                                                 -  BETAS.^2/(2* ...
                                                  sigma_beta^2));

