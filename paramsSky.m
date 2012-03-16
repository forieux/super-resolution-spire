alpha_step = 2;
beta_step = 2;

sky_alpha_period = alpha_step;
sky_beta_period  = sky_alpha_period;

%% 0.6 is the best to approx the equivalent sinc function...
sigma_alpha = 0.6*sky_alpha_period;
sigma_beta  = 0.6*sky_beta_period;

%% The number of order (so 1 if you want only order 0)
Norder = 1;

%% Gaussian decomposition function
alphaS = [-5:5];
[ALPHAS BETAS] = ndgrid(alphaS,alphaS);

fgaussian = 1/(2*pi*sigma_alpha*sigma_beta)* ...
    exp(-ALPHAS.^2/(2*sigma_alpha^2) -  BETAS.^2/(2* sigma_beta^2));

