function grad = calcQuadGrad(object, data, hypers, Hrond, index, ...
                             regOps, objectMean, Nalpha, Nbeta, Norder, ...
                             Nbolo, Nspeed, Nscan, uspeed, ...
                             theSpeeds, varargin)
  %% CALCQUADGRAD - Compute the gradient of the quadratic criterion at
  %% object
  %%
  %% grad = calcQuadGrad(object, data, hypers, bands, Hrond250,
  %% Hrond360, Hrond520, index250, index360, index520, regOps,
  %% objectMean, Nalpha, Nbeta, Norder)
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
  %% Hrond* -- the transfert function for each band in the
  %% directInvariant convention.
  %%
  %% index* -- the index of observed pixel for each band in the
  %% directInvariant convention.
  %%
  %% regOps -- a Nalpha x Nbeta x N tab that contains the N
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
  %% grad = calcQuadGrad(object, data, hypers, bands, Hrond250,
  %% Hrond360, Hrond520, index250, index360, index520, regOps,
  %% objectMean, Nalpha, Nbeta, Norder)
    
  %% Extras options
  Nbaseoption = 15;
  if nargin > Nbaseoption
    for iargin = 1:2:nargin-Nbaseoption
      if strcmp(varargin{iargin},'hes')
        theHessian = varargin{iargin+1};
      elseif strcmp(varargin{iargin},'fig')
        numfig = varargin{iargin+1};
      end
    end
    
  end    
  
  %% Data reproduction
  output = directInvariant(object, Hrond, index, Nalpha, Nbeta, Norder, Nbolo, Nspeed, Nscan, uspeed, theSpeeds);
  
  %% Compute the diff between model reproduction and data
  error = cell(size(data));
  
  %% 'e = y - Hx'
  for iscan = 1:Nscan
    error(iscan) = {output{iscan} - data{iscan}};
  end
  
  %% Transpose the diff 'H^t e'
  adeq = transposeInvariant(error, Hrond, index, Nalpha, Nbeta, Norder, Nbolo, Nspeed, Nscan, uspeed, theSpeeds);
  
  %% '\gammaB*H^t e'
  adeq = hypers(1,1)*adeq;
    
  %% Regularization
    
  %% Compte Qx in Fourier space (F^\dag Lambda_Q Fx) object(:,:,1:3:end)
  %% is all the order of 250 (take each 3)
  regGrad = calcRegPartGrad(object, objectMean, ...
                            squeeze(hypers(2:end,:)), regOps, Norder);

  %% Full gradient
  grad = 2*adeq + 2*regGrad;
  
  %% Correction by the inverse of the hessian He. So the gradient
  %% direction become He^-1 d_k = F^dag Lambda_He^-1 F d_k
  if exist('theHessian','var')
    grad = real(uifft2(ufft2(grad)./theHessian));
  end
  
  %% Plotting
  if exist('numfig','var')
    
    sfigure(numfig);
    
    subplot(2,2,1)
    imagesc(object);
    axis image
    axis xy
    colormap(gray)

    title('Image')

    subplot(2,2,2)
    imagesc(grad);
    axis image
    axis xy
    colormap(gray)

    title('\Delta', 'interpreter', 'tex')
        
    subplot(2,2,3)
    plot(object(200,:))
    title('Image(0)(\beta)', 'interpreter', 'tex')
    subplot(2,2,4)
    plot(object(:,350))
    title('Image(0)(\alpha)', 'interpreter', 'tex')
                
    drawnow
    
   end

end

function grad = calcRegPartGrad(object, objectMean, hypers, regOps, ...
                                Norder)
  %% CALCREGPARTGRAD - Compute the gradient of the quadratic
  %% regularization criterion
  %%
  %% grad = calcRegPartGrad(object, hypers, objectMean, regOps,  Norder)
  %%
  %% compute the regularization part gradient of the criterion of object
  %% for one band. This function is indep of the operator. A good choice
  %% is the penalization on the continus function. The correlation can
  %% be computed in Fourier space with heuristique approach. It returns
  %%
  %% reg = 2 gamma D^t D x = 2 gamma Q x = 2 gamma F^dag Lambda_Q Fx
  %%
  %% The input parameter is Lambda_Q.
  %%
  %% The parameters regOps can be an transfert function in fourrier
  %% space or a kernel in direct space. The code decide it is a kernel
  %% if the dim is different from the image. It decide it is a transfert
  %% function, and do the operation in fourrier space if the dim are the
  %% same (Nalpha x Nbeta).
  %%
  %% INPUT PARAMETERS
  %%
  %% object -- the object must be an Nalpha x Nbeta x Norder. An object
  %% for one band
  %%
  %% hypers -- the hyperparamter value tab of N x Norder with the number
  %% of line is the number of regularization operator and Norder column.
  %%
  %% objectMean -- the mean of the object law. Must be an Nalpha x Nbeta
  %% x Norder
  %%
  %% regOps -- a Nalpha x Nbeta x N x Norder tab that contains the N
  %% regularization operators (a diff operator, a mean operator etc...)
  %% for each order in fourrier space. OR a supportAlpha x supportBeta x
  %% N x Norder regularization kernel in direct space.
  %%
  %% Norder -- the number of decomposition of lambda (Norder == 1 imply
  %% there is only order 0);
  %%
  %% FUNCTION CALL
  %%
  %% grad = calcRegPartGrad(object, hypers, objectMean, regOps,  Norder)
   
  grad = zeros(size(object));
  
  if size(regOps,1) == size(object,1)
    
    for iorder = 1:Norder
      xrond = ufft2(object(:,:,iorder) - objectMean(:,:,iorder));
      
      for ioperator = 1:size(regOps,3)
        %% Simply F^dag Lambda_Q Fx
        contrib = real(uifft2(regOps(:,:,ioperator,iorder).*xrond));
        grad(:,:,iorder) = grad(:,:,iorder) + hypers(ioperator, ...
                                                     iorder)*contrib;
      end
      
    end
    
  else
    
    for iorder = 1:Norder
      
      im = object(:,:,iorder) - objectMean(:,:,iorder);
      
      for ioperator = 1:size(regOps,3)
        %% Simply Q x which is a convolution
        contrib = conv2(im,regOps(:,:,ioperator,iorder),'same');
        grad(:,:,iorder) = grad(:,:,iorder) + hypers(ioperator, ...
                                                     iorder)*contrib;
      end
      
    end
    
  end
  
end
