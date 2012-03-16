function diffOp = diffOpAlpha(sigmaAlpha,sigmaBeta,stepAlpha, stepBeta, N)
%% DIFFOPALPHA - Compute difference operator for function df(a,b)/db
%%
%% Compute the difference opertor on coefficient in alpha direction for one
%% order when the function is decomposed on gaussian.
%% 
%% It returns the first line of the correlation matrix Q in alpha in x^tQx
%% where f(x) is decomposed on gaussian with coefficients x.
%% 
%% FUNCTION CALL
%% 
%% diffOp = diffOpAlpha(sigmaAlpha,sigmaBeta,stepAlpha,stepBeta)
%% 
%% PARAMETERS
%% 
%% sigmaAlpha :  standart deviation of the gaussian in alpha
%% sigmaBeta  :  standart deviation of the gaussian in beta
%% stepAlpha  :  step between 2 succesive gaussian in alpha
%% stepBeta   :  step between 2 succesive gaussian in beta
%% N          :  number of sigma (N*sigma) for the support size

%% Total sigma of the gaussian
sigmaTotA = 2*sigmaAlpha^2 / stepAlpha^2;
sigmaTotB = 2*sigmaBeta^2 / stepBeta^2;

%% Support of the kernel
alpha = -round(N*sqrt(sigmaTotA)) : round(N*sqrt(sigmaTotA));
beta = -round(N*sqrt(sigmaTotB)) : round(N*sqrt(sigmaTotB));
[ALPHAS,BETAS] = ndgrid(alpha,beta);

%% The kernel
normalization = 8*pi*sigmaAlpha^5*sigmaBeta;
kernel = ( sigmaAlpha^2 - stepAlpha^2 / 2 * ALPHAS.^2 ) .* exp( -( ...
    ALPHAS.^2/sigmaTotA + BETAS.^2/sigmaTotB ) / 2 ) / normalization;

diffOp = kernel;

end



