function grad = appHessian2(object, hypers, Hdirect, coefs, diffOps, ...
                           Nalpha, Nbeta, Norder, Nspeed)
  %% CALCQUADGRAD - Compute the gradient of the quadratic criterion at
  %% object
  %%
  %% grad = appHessian(object, data, hypers, bands, Hdirect250, Hdirect360,
  %% Hdirect520, coefs250, coefs360, coefs520, diffOps, objectMean, Nalpha,
  %%  Nbeta, Norder)
  %%
  %% compute the criterion gradiant at object. Data adequation and
  %% regularization are quadratic. This function is indep of the
  %% operator. A good choice is to make the penalization on the continus
  %% function.  The correlation is computed in Fourier space. It returns
  %%
  %% grad = 2 gammaB H^t(y - Hx) + 2 gamma D^t D x = 2 gammaB H^t(y -
  %% Hx) + 2 gamma Q x = 2 gammaB H^t(y - Hx) + 2 gamma F^dag Lambda_Q
  %% Fx
  %%
  %% The input parameter is Lambda_Q and necessery parameter to compute
  %% Hx and H^t y
  %%
  %% calcQuadGrad(..., 'fig', NUMFIG) plot in addition some element on
  %% figure NUMFIG
  %%
  %% calcQuadGrad(..., 'hes', HES) make a correction of the gradient of
  %% order 0 by the inverse of the hessian HES of order 0 in fourrier
  %% space. HES must be the hessian of order 0 in Fourier space.
  %%
  %% INPUT PARAMETERS
  %%
  %% object -- the object must be an Nalpha x Nbeta x (3*Norder) tab
  %% ordered in 250, 360 and 520.
  %%
  %% data -- are the data in directInvariant output convention. A three
  %% cell one for each band. Each cell contains data for all the scan
  %% one cell for each scan.
  %%
  %% hypers -- the hyperparamter value tab of (N+1) x 3 x Norder with
  %% the number of line is the number of regularization operator, for
  %% each band in column, and the third dimension for the order. The
  %% first line hypers(1,:,1) is the noise precision (so indep of order)
  %% for each band.
  %%
  %% bands -- a vector of 3 with 1 to compute a band ordered in 250, 360
  %% and 520. Ex : [0 1 0] compute only for 360.
  %%
  %% Hdirect* -- the transfert function for each band in the
  %% directInvariant convention.
  %%
  %% coefs* -- the index of observed pixel for each band in the
  %% directInvariant convention.
  %%
  %% diffOps -- a Nalpha x Nbeta x N tab that contains the N
  %% regularization operators (a diff operator, a mean operator etc...).
  %% The same is apply at each order
  %%
  %% objectMean -- the mean of the object law.  Must be an Nalpha x
  %% Nbeta x (3*Norder) tab ordered in 250, 360 then 520.
  %%
  %% Nalpha, Nbeta, Norder -- the number of alpha, beta and
  %% decomposition of lambda (Norder == 1 imply there is only order 0);
  %%
  %% FUNCTION CALL
  %%
  %% grad = calcQuadGrad(object, data, hypers, bands, Hdirect250,
  %% Hdirect360, Hdirect520, index250, index360, index520, diffOps,
  %% objectMean, Nalpha, Nbeta, Norder)
  
  x_tmp = zeros(Nalpha, Nbeta, Nspeed);
  dataAdeq = zeros(Nalpha, Nbeta, Norder);
  priorAdeq = zeros(size(dataAdeq));

  %% data adequation
  %% Make a convolution for each order for each speed
  for iorder = 0:Norder-1

    for ispeed = 1:Nspeed
      
      %% The convolution
      x_iorder = conv2(object(:,:,iorder*3+1), ...
                       Hdirect(:,:,iorder+1,ispeed), 'same');
      
      %% Sum of order (but depend on speed)
      x_tmp(:,:,ispeed) = x_tmp(:,:,ispeed) + x_iorder;
      
    end
    
  end

  %% Decimation
  x_tmp = x_tmp.*coefs;

  %% For each speed and each order, to the convolution in Fourier space
  %% with the transpose or the conjugate.
  for ispeed = 1:Nspeed
    
    for iorder = 0:Norder-1

      %% Convolution in Fourier space and addition for one order off all
      %% scan even with different TF.
      obj = conv2(x_tmp(:,:,ispeed), fliplr(flipud(Hdirect(:,:, iorder+1, ispeed))), 'same');
      dataAdeq(:,:,iorder*3+1) = dataAdeq(:,:,iorder*3+1) + obj;

    end    
    
  end
  
  for iorder = 1:Norder
    dataAdeq(:,:,iorder) = hypers(1,iorder)*dataAdeq(:,:,iorder);
  end

  %% Regularization
    
  %% Compte Qx in Fourier space (F^\dag Lambda_Q Fx) object(:,:,1:3:end)
  %% is all the order of 250 (take each 3)
  for iorder = 1:Norder
    
    for ioperator = 1:size(diffOps,3)
      contrib = conv2(object(:,:,iorder), diffOps(:,:,ioperator,iorder), 'same');
      priorAdeq(:,:,iorder) = priorAdeq(:,:,iorder) + ...
          hypers(ioperator+1, iorder)*contrib;
    end
    
  end
    
  %%

  %% Full gradient
  grad = 2*(dataAdeq + priorAdeq);
    
end

