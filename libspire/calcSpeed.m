function [uspeed, theSpeeds, Nspeed] = calcSpeed(pointing, Nscan, Nbolo, tPeriod)

  theSpeeds = zeros(2,1,Nscan);
  
  %% Projection
  for iscan = 1:Nscan
    pScan = pointing{iscan};
    
    pAlpha = reshape(pScan(1,:), size(pScan(1,:),2)/Nbolo, Nbolo);
    pBeta = reshape(pScan(2,:), size(pScan(2,:),2)/Nbolo, Nbolo);
    
    diffAlpha = pAlpha(2:end,:) - pAlpha(1:end-1,:);
    diffBeta = pBeta(2:end,:) - pBeta(1:end-1,:);
    
    theSpeeds(1,1,iscan) = mean(mean(diffAlpha(:)))/tPeriod;
    theSpeeds(2,1,iscan)= mean(mean(diffBeta(:)))/tPeriod;
  end

  speeds_sign = sign(theSpeeds);
    
  uspeed = zeros(2,4);
  uspeed_sign = [1 -1 1 -1; 1 -1 -1 1];

  uspeed(1,1:2) = mean(abs(theSpeeds(1,1,1:Nscan/2)));
  uspeed(2,1:2) = mean(abs(theSpeeds(2,1,1:Nscan/2)));
  uspeed(1,3:4) = mean(abs(theSpeeds(1,1,Nscan/2+1:Nscan)));
  uspeed(2,3:4) = mean(abs(theSpeeds(2,1,Nscan/2+1:Nscan)));
  
  uspeed = uspeed.*uspeed_sign;
  for idx = 1:Nscan
    temp = speeds_sign(:,:,idx)'*uspeed_sign;
    ind = find(temp == 2);
    theSpeeds(:,1,idx) = uspeed(:,ind);
  end
  
  Nspeed = length(uspeed(1,:))

end
