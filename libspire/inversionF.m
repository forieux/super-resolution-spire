function [map ooptim] = inversionF(init, cgoptions, data, hypers, ...
                                   Hrond, index, coefs, offsets, ...
                                   regOp, objectMean, Nalpha, Nbeta, ...
                                   Norder, Nscan, Nbolo, Nspeed, ...
                                   uspeed, theSpeeds)
  %% INVERSION2
  %%
  %%

  %% Sum of coefficients per speed for application of Hessian
  scoefs = zeros(Nalpha, Nbeta, Nspeed);
  for iscan = 1:Nscan
    sIndex = find(uspeed(1,:) == theSpeeds(1,1,iscan));
    scoefs(:,:,sIndex) = scoefs(:,:,sIndex) + coefs{iscan};
  end
  
  %% Correction of offset
  ndata = cell(size(data));
  for iscan = 1:Nscan
    ndata(iscan) = {data{iscan} - repmat(offsets, ...
                                         size(data{iscan},1),1)};
  end
  
  %% Retro-projection of data in Fourier sky space
  dataProj = transposeInvariantF(ndata, Hrond, index, Nalpha, Nbeta, ...
                                 Norder, Nbolo, Nspeed, Nscan, ...
                                 uspeed, theSpeeds);
  for iorder = 1:Norder
    dataProj(:,:,iorder) = 2*hypers(1,iorder)*dataProj(:,:,iorder);
  end

  %% 2*Q*m comming from gradient of (x-m)^tQ(x-m)
  for iRegOp = 1:size(regOp,3)
    for iorder = 1:Norder
      dataProj = dataProj + 2*hypers(iRegOp+1,iorder)*regOp(:,:,iRegOp).*objectMean;
    end
  end
  
  %% Optimisation
  [map sortieops histo] = conjGrad('appHessian', init, dataProj, ...
                                   cgoptions, hypers, Hrond, scoefs, ...
                                   regOp, Nalpha, Nbeta, Norder, ...
                                   Nspeed);

  ooptim = {sortieops histo};
  
end
