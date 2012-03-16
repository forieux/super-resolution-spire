function crit = calcQuadCrit(object, data, hypers, Hrond, index, ...
                             regOps, objectMean, Nalpha, Nbeta, ...
                             Norder, Nbolo, Nspeed, Nscan, uspeed, ...
                             theSpeeds, varargin)
  %% CALCQUADCRIT - Compute the value of the quadratic criterion value
  %% of object
  %%
  %% val = calcQuadCrit(object, data, hypers, bands, Hrond250, Hrond360,
  %%  Hrond520, index250, index360, index520, regOps, objectMean,
  %%  Nalpha, Nbeta, Norder)
  %%
  %% compute the criterion value at object. Data adequation and
  %% regularization are quadratic. This function is indep of the
  %% correlation. A good choice is to make the penalization on the
  %% continous function. The correlation is computed in Fourier space.
  %%
  %% INPUT PARAMETERS
  %%
  %% object -- the object must be an Nalpha x Nbeta x (3*Norder) tab
  %% ordered in 250, 360 then 520.
  %%
  %% data -- are the data in directInvariant output convention. A three
  %% cell one for each band. Each cell contains data for all the scan
  %% one cell for each scan.
  %%
  %% hypers -- the hypersparamter value tab of (N+1) x 3 x Norder with
  %% the number of line is the number of regularization operator, for
  %% each band in column, and the third dimension for the order. The
  %% first line hypers(1,:,1) is the noise precision (so indep of order)
  %% for each band.
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
  %% regOps -- a Nalpha x Nbeta x N x Norder tab that contains the N
  %% regularization operators (a diff operator, a mean operator etc...)
  %% for each order in Fourier space.
  %%
  %% objectMean -- the mean of the object law.  Must be an Nalpha x
  %% Nbeta x (3*Norder) tab ordered in 250, 360 then 520.
  %%
  %% Nalpha, Nbeta, Norder -- the number of alpha, beta and
  %% decomposition of lambda (Norder == 1 imply there is only order 0);
  %%
  %% FUNCTION CALL
  %%
  %% crit = calcQuadCrit(object, data, hypers, bands, Hrond250,
  %% Hrond360, Hrond520, index250, index360, index520, regOps,
  %% objectMean, Nalpha, Nbeta, Norder)

  %% Data reproduction Hx
  output = directInvariant(object, Hrond, index, Nalpha, Nbeta, ...
                           Norder, Nbolo, Nspeed, Nscan, uspeed, ...
                           theSpeeds);
  
  adeq = 0;
  
  %% ||y - Hx||^2
  for iscan = 1:Nscan
    adeq = adeq + sum(sum((data{iscan} - output{iscan}).^2));
  end
  
  %% gammaB*||y - Hx||^2
  crit = hypers(1,1)*adeq;
  
  %% Regularization gammaX*x^tDx for each band
  reg = calcRegPart(object, objectMean, squeeze(hypers(2:end,:)), regOps, Norder);
  crit = crit + reg;
    
end

function crit = calcRegPart(object, objectMean, hypers, regOps, Norder)
  %% CALCREGPART - Compute the value of the quadratic regularization
  %% criterion
  %%
  %% val = calcRegPart(object, objectMean, hypers, regOps, Norder)
  %%
  %% compute the regularization part of the criterion of object for one
  %% band. This function is indep of the operator. A good choice is the
  %% penalization on the continus function. The correlation is computed
  %% in Fourier space.
  %%
  %% INPUT PARAMETERS
  %%
  %% object -- the object must be an Nalpha x Nbeta x Norder
  %%
  %% objectMean -- the mean of the object law. Must be an Nalpha x Nbeta
  %% x Norder
  %%
  %% hypers -- the hyperparamter value tab of N x Norder with the number
  %% of line is the number of regularization operator and Norder column.
  %%
  %% regOps -- a Nalpha x Nbeta x N x Norder tab that contains the N
  %% regularization operators (a diff operator, a mean operator etc...)
  %% for each order in Fourier space
  %%
  %% Norder -- the number of decomposition of lambda (Norder == 1 imply
  %% there is only order 0);
  %%
  %% FUNCTION CALL
  %%
  %% val = calcRegPart(object, objectMean, hypers, regOps, Norder)
  
  crit = 0;
  
  if size(regOps,1) == size(object,1)
    %% In fourrier space !
    for iorder = 1:Norder
      %% compute Fx
      xiorder = ufft2(object(:,:,iorder) - objectMean(:,:,iorder));
      
      for ioperator = 1:size(regOps,3)
        %% gammaI (Fx^dag) Q_I (Fx)
        contrib = real(conj(xiorder).*regOps(:,:,ioperator,iorder).* ...
                       xiorder);
        crit = crit + hypers(ioperator,iorder)*sum(contrib(:)); 
      end
      
    end
    
  else
    %% In direct space !        
    for iorder = 1:Norder
      
      im = object(:,:,iorder) - objectMean(:,:,iorder);
      
      for ioperator = 1:size(regOps,3)
        reg = regOps(:,:,ioperator,iorder);
        %% Simply x Q x where c = Qx is a convolution and xc a
                                %terme wise product.
        contrib = im.*conv2(object,reg,'same');
        crit = crit + hypers(ioperator, iorder)*sum(contrib(:));
      end
      
    end
    
  end
  
end

