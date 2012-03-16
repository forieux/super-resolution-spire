function object = wiener0(dataRond, hypers, bands, Hrond250, ...
                          Hrond360, Hrond520, regOps, Nalpha, Nbeta, ...
                          varargin)
  %% WIENER0 - Wiener/Hunt filter for order zero
  %%
  %% This code compute the Wiener/Hunt filter for order zero of the sky.
  %% Use the hessian approximation.
  %%
  %% FUNCTION CALL
  %%
  %% object = wiener0(dataRond, hypers, bands, Hrond250, Hrond360,
  %% Hrond520, regOps, Nalpha, Nbeta)
  %%
  %% PARAMETERS
  %%
  %% dataRond -- the NORMALIZED fourrier transform of retroprojected of
  %% data on sky space, with a mean by the number of time the pixel as
  %% been observed for each. This must be a 3D TAB of Nalpha x Nbeta x 3
  %% one for 250, 360 and 520 in that order. Here Norder is supposed to
  %% be 1 (only order 0).
  %%
  %% hypers -- the hypersparamter value tab of (N+1) x 3 x Norder with
  %% the number of line is the number of regularization operator, for
  %% each band in column, and the third dimension for the order. The
  %% first line hypers(1,:,1) is the noise precision (so indep of order)
  %% for each band.
  %%
  %% bands -- is vector of dim 3. If Params(1) equal to 1, the transpose
  %% for 250 is computed. Params(2) and Params(3) for 360 and 520
  %% respectively.
  %%
  %% Hrond250 (360/520) -- a tab of Nalpha x Nbeta x Norder x Nspeed
  %% that contains the DIRECT (the conjugate is automaticlly use)
  %% transfert function for 250 (360/520). Nalpha and Nbeta are the
  %% number of pixel in alpha and beta, respectively. Norder is the
  %% number of order, Nspeed is the number of speed (typicaly four).
  %%
  %% regOps -- a Nalpha x Nbeta x N x Norder tab that contains the N
  %% regularization operators (a diff operator, a mean operator etc...)
  %% for each order in Fourier space.
  %%
  %% Nalpha, Nbeta -- the number of pixel in alpha and beta.
  %%
  %% FUNCTION CALL
  %%
  %% object = wiener0(dataRond, bands, gammaB, Hrond250, Hrond360,
  %% Hrond520, paramsReg, regOp, Nalpha, Nbeta)

    %% Init
    paramsInstrument
    paramsObservation
    
    compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);
    
    object = zeros(Nalpha, Nbeta, 3);
    
    %% Compute the Hessian
    denom = hessian(hypers, bands, Hrond250, Hrond360, Hrond520, regOps, ...
                    Nalpha, Nbeta, 0, 0)./2;
    
    %% Quite simple in fact. The Hessian is (gammaB H^tH + gammaXD^tD)
    %% in Fourier space. Compute the Wiener/Hunt filter
    
    %% F^\dag gammaB * (gammaB H^tH + gammaXD^tD)^1 conj(H^t) y

    if compute250
        object(:,:,1) = real(uifft2(sum(hypers(1,1,1)*conj(Hrond250(:,:, ...
                                                          1,:)),4).* ...
                                     dataRond(:,:,1)./denom(:,:,1)));
    end
    
    if compute360
        object(:,:,2) = real(uifft2(sum(hypers(1,2,1)*conj(Hrond360(:,:, ...
                                                          1,:)),4).* ...
                                     dataRond(:,:,2)./denom(:,:,2)));
    end
    
    if compute520
        object(:,:,3) = real(uifft2(sum(hypers(1,3,1)*conj(Hrond520(:,:, ...
                                                          1,:)),4).* ...
                                     dataRond(:,:,3)./denom(:,:,3)));
    end 
    
end
