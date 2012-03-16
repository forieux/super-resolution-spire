function sky = createSky(file, Nalpha, Nbeta, Norder,flux)
%% SIMULATEDATA - Simulate data
    
    load(file);
    %% to adapte the size
    %% Repmat the sky for each band
    sky = repmat(sky, [1 1 3*Norder]);
    %% Futur sky with this dim
    sky2 = zeros(Nalpha, Nbeta, 3*Norder);
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
     
end
