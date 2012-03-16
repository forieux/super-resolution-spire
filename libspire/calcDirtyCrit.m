function crit = calcDirtyCrit(object, data, hypers, bands, index250, ...
                              index360, index520, Nalpha, Nbeta)
  %% CALCQUADCRIT - Compute the value of the quadratic criterion value
  %% of object
  %%
  %% val = calcQuadCrit(object, data, hypers, bands, index250, index360,
  %% index520, regOps, Nalpha, Nbeta)
  %%
  %% INPUT PARAMETERS
  %%
  %% object -- the object must be an Nalpha x Nbeta x (3*Norder) tab
  %% ordered in 250, 360 and 520.
  %%
  %% data -- are the data in directDirty output convention. A three cell
  %% one for each band. Each cell contains data for all the scan one
  %% cell for each scan.
  %%
  %% hypers -- the hypersparamter value tab of 2 x 3 for each band in
  %% column. The first line hypers(1,:,1) is the noise precision. The
  %% second line is for diff.
  %%
  %% bands -- a vector of 3 with 1 to compute a band ordered in 250, 360
  %% and 520. Ex : [0 1 0] compute only for 360.
  %%
  %% index* -- the index of observed pixel for each band in the
  %% directInvariant convention.
  %%
  %% Nalpha, Nbeta -- the number of alpha, beta
  %%
  %% FUNCTION CALL
  %%
  %% val = calcDirtyCrit(object, hypers, data, bands, index250,
  %% index360, index520, Nalpha, Nbeta)

  %% Init
  paramsInstrument
  paramsObservation
  
  compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);    
  
  %% Data adequation
  output = directDirty(object, bands, index250, index360, index520, ...
                       Nalpha, Nbeta);
  
  data250 = data{1}; data360 = data{2}; data520 = data{3};
  sortie250 = output{1}; sortie360 = output{2}; sortie520 = output{3};
  
  adeq250 = 0; adeq360 = 0; adeq520 = 0;
  for iscan = 1:N_scan_total
    if compute250
      adeq250 = adeq250 + sum(sum((data250{iscan} - sortie250{iscan}).^2));
    end
    if compute360
      adeq360 = adeq360 + sum(sum((data360{iscan} - sortie360{iscan}).^2));
    end
    if compute520
      adeq520 = adeq520 + sum(sum((data520{iscan} - sortie520{iscan}).^2));
    end
  end
  
  crit = hypers(1,1)*adeq250 + hypers(1,2)*adeq360 + hypers(1,3)* adeq520;
  
  Doux250 = 0; Doux360 = 0; Doux520 = 0;
  %% Regularization
  if compute250
    DouxAlpha = diff(object(:,:,1),1,1).^2;
    DouxBeta = diff(object(:,:,1),1,2).^2;
    Doux250 = sum(DouxAlpha(:)) + sum(DouxBeta(:));
  end    
  if compute360
    DouxAlpha = diff(object(:,:,2),1,1).^2;
    DouxBeta = diff(object(:,:,2),1,2).^2;
    Doux360 = sum(DouxAlpha(:)) + sum(DouxBeta(:));
  end    
  if compute520
    DouxAlpha = diff(object(:,:,3),1,1).^2;
    DouxBeta = diff(object(:,:,3),1,2).^2;
    Doux520 = sum(DouxAlpha(:)) + sum(DouxBeta(:));
  end    
  
  crit = crit + hypers(2,1)*Doux250 + hypers(2,2)*Doux360 + hypers(2,3)*Doux520;
  
end
