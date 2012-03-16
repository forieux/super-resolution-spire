function [minimizer varargout] = pconjGrad(fcalcAx, init, dataProjected, ...
                                options, fpreCond, varargin)
  %% PRECONDCONJGRAD - Preconditionner Conjugate gradient optimisation
  %%
  %% [m] = conjGrad('calcAx', init, b, options, ...)
  %%
  %% provide a minimizer 'm' of the criterion J(x) = 1/2 x^tAx + b^tx as
  %% the solution of the linear system Ax = b, computed by a conjugate
  %% gradient descente algorithm.  This implementation is adapted, but
  %% absolutely not restricted, to inverse problems where the criterion
  %% take the form
  %%
  %% J(x) = ||y - Hx||^2 + l||Dx||^2, and a algorithm is available to
  %% compute Hx or H^te (x is the unkown, H the direct model, y is data
  %% and D the regularisation.
  %%
  %% In these cases, A = 2(H^tH + lD^tD) and b = H^ty.
  %%
  %% Anyway, the code is quite general since it accept any function that
  %% compute the product Ax, whatever A.
  %%
  %% PARAMETERS
  %%
  %% 'calcAx' - the name of the function to compute the matrix vector
  %% product Ax. Consequently, it's not necessary to compute and store
  %% matrix A, only the product is necessary (think about the
  %% convolution). The name can be anything, but must correspond to a
  %% callable function. The first argument of 'calcAx' MUST BE the
  %% vector x.
  %%
  %% init - the starting point of the optimisation x^(0)
  %%
  %% b - the term b = H^t y. Must be provided.
  %%
  %% options - the options of the algorithm. See section below.
  %%
  %% ... - all the remaning argument are passed to the function 'calcAx'
  %% which is call like this calcAx(x,...)
  %%
  %% OPTIONS
  %%
  %% The variable options is matlab structure. The field are listed
  %% below. Options with '*' are necessary. Others are optional.
  %%
  %% thresold (*) - the stoping criterion
  %%
  %% maxIter (*) - the maximum number of iteration (the algorithm is
  %% automatically stopped when the iteration equal the dimension of x)
  %%
  %% numfig - this option must be a integer. In these case, the current
  %% minimizer and residual are displayed as image of real number in the
  %% figure 'numfig'.
  %%
  %% An example is
  %%
  %% cgoptions.thresold = 1e-6; cgoptions.maxIter = 50;

  try numfig = options.numfig;, end;
  try fcalcCrit = options.fcalcCrit;, end;
  try fcalcCritArgs = options.fcalcCritArgs;, end;

  %% Gradient at current init
  residual = dataProjected - feval(fcalcAx, init, varargin{:});
  descente = feval(fpreCond,residual, varargin{:});
  resNormInit = real(sum(conj(residual(:).*descente(:))));
  secant = descente;
  
  minimizer = init;
  histo = zeros(1,3);
  
  %% Compute the criterion value
  if exist('fcalcCrit','var')
    if exist('fcalcCritArgs','var')
      histo(1,3) = feval(fcalcCrit,init,varargin{:},fcalcCritArgs{:});
    else
      histo(1,3) = feval(fcalcCrit,init,varargin{:});
    end
  end
  
  for iIter = 2:numel(init)

    T0 = cputime; % For computation of cputime
        
    %% Ad
    q = feval(fcalcAx, descente, varargin{:});
    
    %% a = r^tr/d^tAd
    %% Optimal step in direction of descente
    step = sum(conj(residual(:)).*secant(:))/sum(conj(descente(:)).*q(:));
    
    %% Descente x^(i+1) = x^(i) + ad
    minimizer = minimizer + step*descente;
    %% r^(i+1) = r^(i) - a*Ad (think residual as gradient in data space)
    newResidual = residual - step*q;
    
    %% Conjugate direction with preconditionner
    newSecant = feval(fpreCond,newResidual, varargin{:});
    betaCoef = sum(conj(newResidual(:)).*newSecant(:))/sum(conj(residual(:)).*secant(:));
    descente = newSecant + betaCoef*descente;

    residual = newResidual;
    secant = newSecant;
    
    if exist('numfig','var')
      sfigure(numfig);
      
      subplot(121)
      imagesc(minimizer)
      colormap(gray); colorbar; axis image
      title(['Current point (iter = ',num2str(iIter),')'])
      
      subplot(122)
      imagesc(residual)
      colormap(gray); colorbar; axis image
      title('Residual')
      
      drawnow
    end

    %% Time to compute the loop
    histo(iIter,2) = cputime - T0;
    
    %% Compute the criterion value
    if exist('fcalcCrit','var')
      if exist('fcalcCritArgs','var')
        histo(iIter,3) = feval(fcalcCrit,minimizer,varargin{:},fcalcCritArgs{:});
      else
        histo(iIter,3) = feval(fcalcCrit,minimizer,varargin{:});
      end
    end

    %% Stopping criterion
%    histo(iIter,1) = sum(abs((step*descente(:))).^2)/sum(abs(minimizer(:)).^2);
    histo(iIter,1) = real(sum(conj(newResidual(:)).*secant(:)));
    if( histo(iIter,1) < options.thresold*resNormInit)
      histo(:,2) = cumsum(histo(:,2));
      state = 0; %% end by criterion
      if nargout == 2
        varargout = {state};
      elseif nargout == 3
        varargout = {state, histo};
      end
      break;
    end
    
    if(iIter > options.maxIter)
      histo(:,2) = cumsum(histo(:,2));
      state = 1; %% end by iteration
      if nargout == 2
        varargout = {state};
      elseif nargout == 3
        varargout = {state, histo};
      end
      break;
    end

  end
end
