function object = wiener01(dataRond, bands, gammaB, Hrond250, ...
                           Hrond360, Hrond520, regParams, regOp, ...
                           Nalpha, Nbeta, varargin)
  %% WIENER01 - Wiener/Hunt filter for order 0 and 1
  %%
  %% Compute the Wiener/Hunt filter for order zero and 1. Use the
  %% hessian approximation.
  %%
  %% FUNCTION CALL
  %%
  %% object = wiener01(dataRond, bands, gammaB, Hrond250, Hrond360,
  %% Hrond520, regParams, regOp, index250, index360, index520, Nalpha,
  %% Nbeta)
  %%
  %% PARAMETERS
  %%
  %% dataRond : the NORMALIZED fourrier transform of retroprojected of
  %% data on sky space, with a mean by the number of time the pixel as
  %% been observed for each. This must be a 3D TAB of Nalpha x Nbeta x 3
  %% one for 250, 360 and 520 in that order. Here Norder is supposed to
  %% be 1 (only order 0).
  %%
  %% bands : is vector of dim 3. If Params(1) equal to 1, the transpose
  %% for 250 is computed. Params(2) and Params(3) for 360 and 520
  %% respectively.
  %%
  %% gammaB, gamma1, gamma0 : the regularization parameters. gammaB is
  %% the inverse variance of the noise, gamma1 for the diff and gamma0
  %% for the mean.
  %%
  %% Hrond250 (360/520) : a tab of Nalpha, Nbeta, Nspeed that contains
  %% the DIRECT (the conjugate is automaticlly use) transfert function
  %% for 250 (360/520). Nalpha and Nbeta are the number of pixel in
  %% alpha and beta, respectively. Nspeed is the number of speed
  %% (typicaly four).
  %%
  %% regParams : a Nx3 for regularisation parameters. One column for
  %% 250, 360 and 520 respectively. Each line correspond to a
  %% regulrization operator.
  %%
  %% regOp : a Nalpha x Nbeta x N tab that contains the N regularization
  %% operators (a diff operator, a mean operator etc...)
  %%
  %% index250 (360/520) : are the index corresponding to the position in
  %% (alpha, beta) when the data as been aquired for 250 (350/520). This
  %% must be index so use potentialy the matlab function sub2ind.
  %%
  %% Nalpha, Nbeta : the number of pixel in alpha ans beta.
  %%
  %% FUNCTION CALL
  %%
  %% object = wiener01(dataRond, bands, gammaB, Hrond250, Hrond360,
  %% Hrond520, paramsReg, regOp, index250, index360, index520, Nalpha,
  %% Nbeta)

  %% Init
    paramsInstrument
    paramsObservation
    
    compute250 = bands(1);
    compute360 = bands(2);
    compute520 = bands(3);
    
    %% Partition of the hessien
    h00 = hessian00(bands, gammaB, Hrond250, Hrond360, Hrond520, ...
                    regParams, regOp, Nalpha, Nbeta);
    h10 = hessian10(bands, gammaB, Hrond250, Hrond360, Hrond520, ...
                    Nalpha, Nbeta);
    h01 = hessian01(bands, gammaB, Hrond250, Hrond360, Hrond520, ...
                    Nalpha, Nbeta);
    h11 = hessian11(bands, gammaB, Hrond250, Hrond360, Hrond520, ...
                    regParams, regOp, Nalpha, Nbeta);
    
    %% Inversion of the hessian (everything is diagonal)
    p00 = 1./(h00 - h01.*h10./h11);
    p11 = 1./(h11 - h10.*h01./h00);
    p01 = -h01.*p11./h00;
    p10 = -h10.*p00./h11;
    
    object = zeros(Nalpha, Nbeta, 6);
    object(:,:,1) = dataRond(:,:,1);
    object(:,:,4) = dataRond(:,:,1);
    object(:,:,2) = dataRond(:,:,2);
    object(:,:,5) = dataRond(:,:,2);
    object(:,:,3) = dataRond(:,:,3);
    object(:,:,6) = dataRond(:,:,3);
    
    if compute250
        object(:,:,1) = gammaB(1)*sum(conj(Hrond250(:,:,1,:)),4).* ...
            object(:,:,1);
        object(:,:,4) = gammaB(1)*sum(conj(Hrond250(:,:,2,:)),4).* ...
            object(:,:,4);
        
        object(:,:,1) = real(uifft2(object(:,:,1).*p00(:,:,1) + ...
                                     object(:,:,4).*p01(:,:,1)));
        object(:,:,4) = real(uifft2(object(:,:,1).*p10(:,:,1) + ...
                                     object(:,:,4).*p11(:,:,1)));
    end
    
    if compute360
        object(:,:,2) = gammaB(2)*sum(conj(Hrond360(:,:,1,:)),4).* ...
            object(:,:,2);
        object(:,:,5) = gammaB(2)*sum(conj(Hrond360(:,:,2,:)),4).* ...
            object(:,:,5);
        
        object(:,:,2) = real(uifft2(object(:,:,2).*p00(:,:,2) + ...
                                     object(:,:,5).*p01(:,:,2)));
        object(:,:,5) = real(uifft2(object(:,:,2).*p10(:,:,2) + ...
                                     object(:,:,5).*p11(:,:,2)));
    end
    
    if compute520
        object(:,:,3) = gammaB(3)*sum(conj(Hrond520(:,:,1,:)),4).* ...
            object(:,:,3);
        object(:,:,6) = gammaB(3)*sum(conj(Hrond520(:,:,2,:)),4).* ...
            object(:,:,6);
        
        object(:,:,3) = real(uifft2(object(:,:,3).*p00(:,:,3) + ...
                                     object(:,:,6).*p01(:,:,3)));
        object(:,:,6) = real(uifft2(object(:,:,3).*p10(:,:,3) + ...
                                     object(:,:,6).*p11(:,:,3)));
    end 
end
