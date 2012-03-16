function grad = calcDirtyGrad(object, data, hypers, bands, index250, ...
                              index360, index520, Nalpha, Nbeta, varargin)
  %% CALCDIRTYGRAD - Compute the gradient of the quadratic criterion at
  %% object
  %%
  %% grad = calcDirtyGrad(object, data, hypers, bands, index250,
  %% index360, index520, Nalpha, Nbeta)
  %%
  %% calcQuadGrad(..., 'fig',NUMFIG) plot in addition some element on
  %% figure NUMFIG
  %%
  %% INPUT PARAMETERS
  %%
  %% object -- the object must be an Nalpha x Nbeta x (3*Norder) tab
  %% ordered in 250, 360 and 520.
  %%
  %% data -- are the data in directInvariant output convention. A three
  %% cell one for each band. Each cell contains data for all the scan
  %% one cell for each scan.
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
  %% grad = calcDirtyGrad(object, data, hypers, bands, index250,
  %% index360, index520, nalpha, nbeta)

  paramsInstrument
  paramsObservation
  
  %% Extras options
  
  if nargin > 10
    numfig = varargin{2};
  end    
  compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);
  
  %% data adequation
  output = directDirty(object, bands, index250, index360, index520, ...
                       Nalpha, Nbeta);
  
  output250 = output{1}; output360 = output{2}; output520 = output{3};
  
  data250 = data{1}; data360 = data{2}; data520 = data{3};
  
  %% Compute the diff between model reproduction and data
  error250 = cell(size(data250));
  error360 = cell(size(data360));
  error520 = cell(size(data520));
  
  for iscan = 1:N_scan_total
    if compute250
      error250(iscan) = {output250{iscan} - data250{iscan}};
    end
    if compute360
      error360(iscan) = {output360{iscan} - data360{iscan}};
    end
    if compute520
      error520(iscan) = {output520{iscan} - data520{iscan}};
    end
  end
  
  %% Reformat The Diff for transpose function
  error = {error250, error360, error520};
  
  %% Transpose the diff
  adeq = transposeDirty(error, bands, index250, index360, index520, ...
                          Nalpha, Nbeta);
  
  adeq(:,:,1:3:end) = hypers(1,1)*adeq(:,:,1:3:end);
  adeq(:,:,2:3:end) = hypers(1,2)*adeq(:,:,2:3:end);
  adeq(:,:,3:3:end) = hypers(1,3)*adeq(:,:,3:3:end);
    
  %% Regularization
    
  regGrad = zeros(size(adeq)); 
  
  if compute250
    DouxAlpha = diff(object(:,:,1),1,1);
    DouxAlpha = cat(1,DouxAlpha(1,:),DouxAlpha,DouxAlpha(end,:));
    DouxAlpha = -2*diff(DouxAlpha,1,1);
  
    DouxBeta = diff(object(:,:,1),1,2);
    DouxBeta = cat(2,DouxBeta(:,1),DouxBeta,DouxBeta(:,end));
    DouxBeta = -2*diff(DouxBeta,1,2);
        
    regGrad(:,:,1) = hypers(2,1)*(DouxAlpha + DouxBeta);
  end
  
  if compute360
    DouxAlpha = diff(object(:,:,2),1,1);
    DouxAlpha = cat(1,DouxAlpha(1,:),DouxAlpha,DouxAlpha(end,:));
    DouxAlpha = -2*diff(DouxAlpha,1,1);
  
    DouxBeta = diff(object(:,:,2),1,2);
    DouxBeta = cat(2,DouxBeta(:,1),DouxBeta,DouxBeta(:,end));
    DouxBeta = -2*diff(DouxBeta,1,2);

    regGrad(:,:,2) = hypers(2,2)*(DouxAlpha + DouxBeta);
  end
  
  if compute520
    DouxAlpha = diff(object(:,:,3),1,1);
    DouxAlpha = cat(1,DouxAlpha(1,:),DouxAlpha,DouxAlpha(end,:));
    DouxAlpha = -2*diff(DouxAlpha,1,1);
  
    DouxBeta = diff(object(:,:,3),1,2);
    DouxBeta = cat(2,DouxBeta(:,1),DouxBeta,DouxBeta(:,end));
    DouxBeta = -2*diff(DouxBeta,1,2);
        
    regGrad(:,:,3) = hypers(2,3)*(DouxAlpha + DouxBeta);
  end
    
  %% Gradient
  grad = 2*adeq + 2*regGrad;
  
  %% Plotting
  if exist('numfig','var')
    
    sfigure(numfig);
    
    subplot(2,3,1)
    imagesc(object(:,:,1))
    axis image
    axis xy
    colormap(gray)
    title('250(0)')
    
    subplot(2,3,2)
    imagesc(object(:,:,2))
    axis image
    axis xy
    colormap(gray)
    title('360(0)')
    
    subplot(2,3,3)
    imagesc(object(:,:,3))
    axis image
    axis xy
    colormap(gray)
    title('520(0)')
    
    subplot(2,3,4)
    plot(object(200,:,1))
    title('250(0)(\beta)', 'interpreter', 'tex')
    
    subplot(2,3,5)
    plot(object(200,:,2))
    title('360(0)(\beta)', 'interpreter', 'tex')
    
    subplot(2,3,6)
    plot(object(200,:,3))
    title('520(0)(\beta)', 'interpreter', 'tex')
    
    drawnow
    
  end
  
end


