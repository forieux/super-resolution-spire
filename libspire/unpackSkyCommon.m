function [sky common] = unpackSkyCommon(vec, Nalpha, Nbeta, Norder)
%% UNPACKSKYCOMMON - UnPack the vec for gpac in sky and common
    
  numberOfSky = Nalpha*Nbeta*Norder*3;
  
  sky = vec(1:numberOfSky);
  common = vec(numberOfSky+1:end);

  %% The two in one column vector
  sky = reshape(sky, [Nalpha Nbeta Norder*3]);
  common = reshape(common, [], 3);

end
