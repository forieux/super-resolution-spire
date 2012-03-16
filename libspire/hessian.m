function hessian = hessian(hypers, bands, Hrond250, Hrond360, Hrond520, ...
                           regOps, Nalpha, Nbeta, oA, oB, varargin)
%% HESSIAN01 - Hessian TF approximation for any cross order 
%%
%% This code compute the hessian Identity approximation of spire in Fourrier
%% space for any cross order A and B
%% 
%% FUNCTION CALL
%% 
%% hessian = hessian(hypers, bands, Hrond250, Hrond360, Hrond520, regOps,
%% Nalpha, Nbeta, orderA, orderB)
%% 
%% PARAMETERS
%% 
%% hypers -- the hypersparamter value tab of (N+1) x 3 x Norder with the number
%% of line is the number of regularization operator, for each band in column,
%% and the third dimension for the order. The first line hypers(1,:,1) is the
%% noise precision (so indep of order) for each band.
%% 
%% bands -- is vector of dim 3. If Params(1) equal to 1, the transpose for 250
%% is computed. Params(2) and Params(3) for 360 and 520 respectively.
%%
%% Hrond250 (360/520) -- a tab of Nalpha x Nbeta x Norder x Nspeed that
%% contains the DIRECT (the conjugate is automaticlly use) transfert function
%% for 250 (360/520). Nalpha and Nbeta are the number of pixel in alpha and
%% beta, respectively. Norder is the number of order, Nspeed is the number of
%% speed (typicaly four).
%%
%% regOps -- a Nalpha x Nbeta x N x Norder tab that contains the N
%% regularization operators (a diff operator, a mean operator etc...) for
%% each order in Fourier space.
%%    
%% Nalpha, Nbeta -- the number of pixel in alpha and beta.
%% 
%% orderA, orderB -- the two order. Be carefull, orderA/orderB is different
%% from orderB/orderA. For cross order 0&1 it must be orderA = 0, orderB = 1.
%%  
%% FUNCTION CALL
%% 
%% hessian = hessian(bands, hypers, Hrond250, Hrond360, Hrond520, regOps,
%% Nalpha, Nbeta, orderA, orderB)

  compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);
  
  hessian = zeros(Nalpha, Nbeta, 3);
  
  %% Sum for all scan the product of conj(orderA).*orderB
  if compute250
    hessian(:,:,1) = hypers(1,1,1)*sum(conj(Hrond250(:,:,oA+1,:)).* ...
                                       Hrond250(:,:, oB+1,:),4);
  end
  
  if compute360
    hessian(:,:,2) = hypers(1,2,1)*sum(conj(Hrond360(:,:,oA+1,:)).* ...
                                       Hrond360(:,:, oB+1,:),4);
  end
  
  if compute520
    hessian(:,:,3) = hypers(1,3,1)*sum(conj(Hrond520(:,:,oA+1,:)).* ...
                                       Hrond520(:,:, oB+1,:),4);
  end 
  
  if oA == oB
    
    for ioperator = 1:size(regOps,3)
      
      if compute250
        hessian(:,:,1) = hessian(:,:,1) + hypers(ioperator+1,1,oA+1)* ...
            regOps(:,:,ioperator,oA+1);
      end
      
      if compute360
        hessian(:,:,2) = hessian(:,:,2) + hypers(ioperator+1,2,oA+1)* ...
            regOps(:,:,ioperator,oA+1);
      end
      
      if compute520
        hessian(:,:,3) = hessian(:,:,3) + hypers(ioperator+1,3,oA+1)* ...
            regOps(:,:,ioperator,oA+1);
      end
      
    end
    
  end
  
  hessian = 2*hessian;
  
end
    

