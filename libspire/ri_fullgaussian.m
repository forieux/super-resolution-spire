function response = ri_fullgaussian(alpha, beta, lambda, ...
    opt_efficiency, sigma, ...
    time_constante, gain, speed)
%% RI_FULLGAUSSIAN - Compute the full gaussian model response of spire
%%
%% response = ri_fullgaussian(alpha, beta, lambda, sigma,
%% time_constante, gain, speed)
%%
%% compute the acquisition model for alpha and beta.
%%
%% This function can be used to compute the impultionnal response of
%% the instrument if the invariance hypothesis is used or to compute
%% the H matrix. It's only depend on support parameter provided.
%%
%% All the response have the same energy in each lambda. To get 1,
%% remember to multiply by the surface of a pixel.
%%
%% The response correspond to the model 1: - The sky is decomposed on
%% equaly spaced gaussian fonction - The airy disc is approximed with
%% a gaussian - The feedhorns respons are gaussian - The bolometer
%% response is linear - Scan map protocol - The response are
%% normalized to 1 for each lambda. But remember to multiply by the
%% size of a pixel to get effectively 1.
%%
%% FONCTION CALL
%%
%% response = ri_fullgaussian(alpha, beta, lambda, sigma,
%% time_constante, gain, speed)
%%
%% PARAMETERS
%%
%% alpha, beta, lambda -- The coordinate vector on wich the response
%% must be computed in [arcsec arcsec meter].
%%
%% opt_efficient -- the global efficiency for optique. Physicaly it is
%% a gain between 0 and 1.
%%
%% sigma -- The vector [coef, orig_alpha, orig_beta] that define the
%% dependence in lambda. This correspond to equations
%%
%%           sigma_alpha = coef*lambda + orig_alpha and sigma_beta =
%%           coef*lambda + orig_beta
%%
%% with lambda in meter and sigma_alpha/beta in arcsec.
%%
%% time_constante -- The time constant in seconde of the bolometer
%% model.
%%
%% gain -- The gain of the bolometer model in ?
%%
%% speed -- The [v_alpha v_beta] vector speed (neccessary for the
%% normalization of the response) in arcsec/seconde.
%%
%% FONCTION CALL
%%
%% response = ri_fullgaussian(alpha, beta, lambda, sigma,
%% time_constante, gain, speed)

  [ALPHAS BETAS LAMBDAS] = ndgrid(alpha, beta, lambda);
  
  coef_sigma = sigma(1);
  orig_sigma_alpha = sigma(2);
  orig_sigma_beta = sigma(3);
  
  v_alpha = speed(1);
  v_beta = speed(2);
  
  SIGMA_alpha = coef_sigma*LAMBDAS + orig_sigma_alpha;
  SIGMA_beta = coef_sigma*LAMBDAS + orig_sigma_beta;
  
  SIGMA_v = sqrt(SIGMA_alpha.^2*v_beta^2 + SIGMA_beta.^2*v_alpha^2);
  
  %% Exponential part indep of the bolometer. In other word it depends on
  %% time only by a shift, or the support.
  H1 = exp(-ALPHAS.^2./(2*SIGMA_alpha.^2) - BETAS.^2./(2*SIGMA_beta.^2));
  
  %% The erfcx part. This the part of influence of the bolometer.
  if time_constante == 0
    H2 = ones(size(H1));
    norm = opt_efficiency*gain./(2*pi.*SIGMA_alpha.*SIGMA_beta);
  else
    arg = SIGMA_alpha.*SIGMA_beta./(sqrt(2).*time_constante.*SIGMA_v) ...
        - ALPHAS.*SIGMA_beta.*v_alpha./(sqrt(2).*SIGMA_alpha.*SIGMA_v) ...
        - BETAS.*SIGMA_alpha.*v_beta./(sqrt(2).*SIGMA_beta.*SIGMA_v);
    H2 = erfcx(arg);
    norm = opt_efficiency*gain./(2*sqrt(2*pi).*SIGMA_v);
  end
  
  response = norm.*H1.*H2;

end
