function norm = normOp(sigmaAlpha, sigmaBeta, stepAlpha, stepBeta)
%% NORMOP - Compute the norm operator for the continus function f(a,b)
%%
%% Compute the norm operator for gaussian decomposed function.
%% 
%% FUNCTION CALL
%% 
%% operator = normOp(sigmaAlpha,sigmaBeta,stepAlpha,stepBeta)
%% 
%% PARAMETERS 
%% 
%% SIGMAALPHA   :  standart deviation of the gaussian in alpha
%% SIGMABETA    :  standart deviation of the gaussian in beta
%% STEPALPHA    :  step between 2 succesive gaussian in alpha
%% STEPBETA     :  step between 2 succesive gaussian in beta

  %% Full standart deviation
  sigmaTotA = 2*sigmaAlpha^2 / stepAlpha^2;
  sigmaTotB = 2*sigmaBeta^2 / stepBeta^2;
  
  %% Support of the kernel
  alpha = -round(10*sqrt(sigmaTotA)) : round(10*sqrt(sigmaTotA));
  beta = -round(10*sqrt(sigmaTotB)) : round(10*sqrt(sigmaTotB));
  [ALPHAS,BETAS] = ndgrid(alpha,beta);
  
  %% Kernel
  normalization = 4*pi*sigmaAlpha*sigmaBeta;
  norm = exp( -( ALPHAS.^2/sigmaTotA + BETAS.^2/sigmaTotB ) / 2 )/ normalization;
  
end
