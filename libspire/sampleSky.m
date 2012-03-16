function [skySample ooptim] = sampleSky(init, options, data, hypers, ...
                                        Hrond, index, coefs, ...
                                        offsets, regOps, Nalpha, ...
                                        Nbeta, Norder, Nscan, Nbolo, ...
                                        Nspeed, uspeed, theSpeeds, ...
                                        varargin)
  %% SAMPLESKY - Simulate the conditionnal posterior law of sky by
  %% optimization
  %%
  %% [skySample SortieOPS Histo] = sampleSky(init, options, data,
  %% hypers, bands, Hrond250, Hrond360, Hrond520, index250, index360,
  %% index520, regOps, Nalpha, Nbeta, Norder)
  %%
  %% compute a sample of the quadratic conditionnal posterior law of the
  %% sky by optimization. the criterion gradiant at object. Data
  %% adequation and regularization are quadratic. This function is indep
  %% of the operator. A good choice is to make the penalization on the
  %% continus function.
  %%
  %% Additionnal output are the results of the optimisation see 'gpac'
  %% help.
  %%
  %% It use the fonction 'calcQuadCrit', 'calcQuadGrad' and 'gpac' to
  %% compute the optimization
  %%
  %% sampleSky(..., 'fig', NUMFIG) plot in addition some element on
  %% figure NUMFIG
  %%
  %% sampleSky(..., 'hes', HES) make a correction of the gradient of
  %% order 0 by the inverse of the hessian HES of order 0 in fourrier
  %% space. HES must be the hessian of order 0 in Fourier space.
  %%
  %% INPUT PARAMETERS
  %%
  %% init -- the init must be an Nalpha x Nbeta x (3*Norder) tab ordered
  %% in 250, 360 and 520.
  %%
  %% options -- the option vector for GPAC, see it's help.
  %%
  %% data -- are the data in directInvariant output convention. A three
  %% cell one for each band. Each cell contains data for all the scan
  %% one cell for each scan.
  %%
  %% hypers -- the hyperparamter value tab of 2 x 3 x Norder. The first
  %% line hypers(1,:,1) is the noise precision (so indep of order) for
  %% each band. The second line is the hypers value for the object.
  %% Column is the band.
  %%
  %% bands -- a vector of 3 with 1 to compute a band ordered in 250, 360
  %% and 520. Ex : [0 1 0] compute only for 360.
  %%
  %% Hrond* -- the transfert function for each band in the
  %% directInvariant convention.
  %%
  %% index* -- the index of observed pixel for each band in the
  %% directInvariant convention.
  %%
  %% regOps -- a Nalpha x Nbeta x Norder tab that contains the inverse
  %% variance (a diff operator, a mean operator etc...).
  %%
  %% Nalpha, Nbeta, Norder -- the number of alpha, beta and
  %% decomposition of lambda (Norder == 1 imply there is only order 0);
  %%
  %% FUNCTION CALL
  %%
  %% sample = sampleSky(init, options, data, hypers, bands, Hrond250,
  %% Hrond360, Hrond520, index250, index360, index520, regOps, Nalpha,
  %% Nbeta, Norder)
    
  stdN = 1/sqrt(hypers(1,1));
  stdX = 1/sqrt(hypers(2,1));

  %% Data perturbation
  for iscan = 1:Nscan
    data(iscan) = {data{iscan} + randn(size(data{iscan}))*stdN};
  end
  
  %% Sky perturbation 
  objectMean = zeros(size(init));
  
  for iorder = 0:Norder-1
    
    if size(regOps,1) == size(init,1)
      %% In fourrier space
      
      %% Standart deviation is the square root of the variance and in
      %% $x^t Q x$, Q is the inverse variance.  Since the transfert
      %% fonction is real (the impultionnal response is even) we can
      %% take the square root of Lambda where Q = F^dag Lambda F
      standartDeviation = 1./sqrt(regOps(:,:,:,iorder+1));
    else
      %% In direct space
      
      %% We have to compute in Fourrier space so compute the
      %% transfert function
      
      inverseVariance = real(ri2fourier(regOps(:,:,:,iorder+1), ...
                                        Nalpha, Nbeta));
      standartDeviation = 1./(inverseVariance.^(1/2));
    end
    standartDeviation(1) = 1;
    
    whiteComplexeNoise = complexeRandn(Nalpha,Nbeta);
    
    objectMean(:,:,iorder+1) = ...
        stdX*standartDeviation.*whiteComplexeNoise;
            
  end
  
  %% Optimization
  [skySample ooptim] = inversionF(init, options, data, hypers, ...
                                  Hrond, index, coefs, offsets, ...
                                  regOps, objectMean, Nalpha, Nbeta, ...
                                  Norder, Nscan, Nbolo, Nspeed, ...
                                  uspeed, theSpeeds);

end

function n = complexeRandn(Nalpha, Nbeta)
%% COMPLEXERANDN - Return a complexe randn sample with correct normalisation
%% 
%% n = complexeRandn(Nalpha, Nbeta)
%% 
%% return the sum a real white idd gaussian sample of std 1/2, plus a
%% imaginary white idd gaussian sample of std 1/2. The sum is normalized by
%% 1/sqrt(number of pixel).
%% 
    
%     n = (randn(Nalpha, Nbeta)*sqrt(1/2) + i*randn(Nalpha, Nbeta)*sqrt(1/2))./ ...
%         sqrt(Nalpha*Nbeta);
    
    n = ufft2(randn(Nalpha, Nbeta));
    
end
