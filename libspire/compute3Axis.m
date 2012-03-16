function [alpha beta] = computeAxis(pointingPSW, pointingPMW, pointingPLW, salpha, sbeta, bound, Nscan)

  minAlpha = 0; minBeta = 0; maxAlpha = 0; maxBeta = 0;
  
  for iscan = 1:Nscan
    p = pointingPSW{iscan};
    minAlpha = min(minAlpha, min(min(p(1,:,:))));
    minBeta  = min(minBeta,  min(min(p(2,:,:))));
    maxAlpha = max(maxAlpha, max(max(p(1,:,:))));
    maxBeta  = max(maxBeta,  max(max(p(2,:,:))));
    p = pointingPMW{iscan};
    minAlpha = min(minAlpha, min(min(p(1,:,:))));
    minBeta  = min(minBeta,  min(min(p(2,:,:))));
    maxAlpha = max(maxAlpha, max(max(p(1,:,:))));
    maxBeta  = max(maxBeta,  max(max(p(2,:,:))));
    p = pointingPLW{iscan};
    minAlpha = min(minAlpha, min(min(p(1,:,:))));
    minBeta  = min(minBeta,  min(min(p(2,:,:))));
    maxAlpha = max(maxAlpha, max(max(p(1,:,:))));
    maxBeta  = max(maxBeta,  max(max(p(2,:,:))));
  end
  
  alpha = minAlpha - bound : salpha : maxAlpha + bound;
  beta = minBeta - bound : sbeta : maxBeta + bound;

  %% The RI will be computed on all the support of the sky, to avoid the
  %% phase problem in Fourrier space. Secondly to be certain that the RI
  %% will be computed on zero coordinate, we need odd number of element.
  if mod(length(alpha),2) == 0
    alpha = [alpha max(alpha) + salpha];
  end
  if mod(length(beta),2) == 0
    beta = [beta max(beta)+sbeta];
  end
  
end
