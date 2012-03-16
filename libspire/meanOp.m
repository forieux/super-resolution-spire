function meanOp = meanOp(stepAlpha, stepBeta, Nalpha, Nbeta)
%% MEANOP - Compute mean operator for function f(a,b)
%%
%% Compute the mean opertor on coefficient one order when the function is
%% decomposed on gaussian.
%% 
%% It returns the mean operator 1 in x^t1x where f(x) is decomposed on
%% gaussian with coefficients x.
%% 
%% FUNCTION CALL
%% 
%% meanOp = meanOp(stepAlpha,stepBeta, Nalpha, Nbeta)
%% 
%% INPUT 
%% 
%% STEPALPHA    :  step between 2 succesive gaussian in alpha
%% STEPBETA     :  step between 2 succesive gaussian in beta
%% Nalpha       :  the number of alpha
%% Nbeta        :  the number of beta
    
    meanOp = ones(Nalpha, Nbeta) ./ (Nalpha*stepAlpha*Nbeta*stepBeta)^2;
    
end
