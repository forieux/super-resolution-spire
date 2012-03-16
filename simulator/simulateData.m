function [data output sky] = simulateData(file, Nalpha, Nbeta, Norder, std, flux, Hrond, index, Nbolo, Nspeed, Nscan, uspeed, theSpeeds)
%% SIMULATEDATA - Simulate data
    
  load(file);
  %% to adapte the size
  %% Repmat the sky for each band
  sky = repmat(sky, [1 1 Norder]);
  %% Futur sky with this dim
  sky2 = zeros(Nalpha, Nbeta, Norder);
  ma = min([Nalpha, size(sky,1)]);
  mb = min([Nbeta,  size(sky,2)]);
  %% Fill with the file
  sky2(1:ma,1:mb,:) = sky(1:ma,1:mb,:);
  %% Put in center
  sky2 = fftshift(circshift(sky2,[-floor(ma/2) -floor(mb/2)]));
  
  sky = sky2; clear sky2;
  
  %% If two order
  if Norder == 2
    for couche = 4:6
      sky(:,:,couche) = -50*flipud(fliplr(sky(:,:,couche)));
    end      
  end
  
  %% In Watt sr^-1 m^-1 (m^-2 ?).
  sky = sky*flux;
  
  output = directInvariant(sky, Hrond, index, Nalpha, Nbeta, Norder, Nbolo, Nspeed, Nscan, uspeed, theSpeeds);

  %% add noise
  for iscan = 1:Nscan
    data(iscan) = {output{iscan} + std*randn(size(output{iscan}))};
  end
  
end
