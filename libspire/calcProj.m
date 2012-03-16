function projection = calcProj(pointing, Nscan)

  temp_ra=[];
  temp_dec=[];
  
  %% determining ra_mean and dec_mean
  for iscan = 1:Nscan
    p = pointing{iscan};
    temp_ra = [temp_ra deg2rad(p(1,:))];
    temp_dec = [temp_dec deg2rad(p(2,:))];
  end
  
  raz = mean(temp_ra);
  decz = mean(temp_dec);
  
  dir_alpha = zeros(1,iscan);
  dir_beta = zeros(1,iscan);
  
  %% Projection
  for iscan = 1:Nscan
    p = pointing{iscan};
    p(2,:) = deg2rad(p(2,:));
    p(1,:) = deg2rad(p(1,:));

    [xi, eta] = convert_ra_dec_bis(p, raz, decz);
  
    projection(iscan) = {[3600.*(180/pi.*xi); rad2arcsec(eta)]};
  end

end
