function [map offsets ooptim] = inversionOffsets(init, cgoptions, data, ...
                                         hypers, Hdirect, index, ...
                                         coefs, offsets, diffOp, ...
                                         objectMean, Nalpha, Nbeta, ...
                                         Norder, Nscan, Nbolo, ...
                                         Nspeed, uspeed, theSpeeds, repeat)
  %% INVERSIONOFFSETS
  %%
  %%

  %% Sum of coefficients per speed for application of Hessian
  scoefs = zeros(Nalpha, Nbeta, Nspeed);
  for iscan = 1:Nscan
    sIndex = find(uspeed(1,:) == theSpeeds(1,iscan));
    scoefs(:,:,sIndex) = scoefs(:,:,sIndex) + coefs{iscan};
  end

  for iRun = 1:repeat
    %% Correction of offset
    ndata = cell(size(data));
    for iscan = 1:Nscan
      ndata(iscan) = {data{iscan} - repmat(offsets, ...
                                           size(data{iscan},1),1)};
    end
    
    %% Retro-projection of data in Fourier sky space
    dataProj = transposeInvariant2(ndata, Hdirect, index, Nalpha, ...
                                   Nbeta, Norder, Nbolo, Nspeed, ...
                                   Nscan, uspeed, theSpeeds);
                                 
    for iorder = 1:Norder
      dataProj(:,:,iorder) = 2*hypers(1,iorder)*dataProj(:,:,iorder);
    end

    %% 2*Q*m comming from gradient of (x-m)^tQ(x-m)
    for iRegOp = 1:size(diffOp,3)
      for iorder = 1:Norder
        dataProj = dataProj + ...
            2*hypers(iRegOp+1,iorder)*conv2(objectMean, ...
                                            diffOp(:,:,iRegOp), 'same');
      end
    end
    
    %% Optimisation
    [map sortieops histo] = conjGrad('appHessian2', init, ...
                                     dataProj, cgoptions, hypers, ...
                                     Hdirect, scoefs, diffOp, ...
                                     Nalpha, Nbeta, Norder, Nspeed);
    
    %% Offsets
    dataRepro = directInvariant2(map, Hdirect, index, Nalpha, Nbeta, ...
                                 Norder, Nbolo, Nspeed, Nscan, ...
                                 uspeed, theSpeeds);
    
    offsets = estimOffsets(data, dataRepro, Nscan);
    
  end
  
  ooptim = {sortieops histo};
  
end
