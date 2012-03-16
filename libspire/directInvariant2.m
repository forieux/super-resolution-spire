function data  = directInvariant2(sky, Hdirect, index, Nalpha, Nbeta, ...
                                 Norder, Nbolo, Nspeed, Nscan, ...
                                 uspeed, theSpeeds)
  %% DIRECTINVARIANT - direct computation of the invariant spire model
  %%
  %% This code compute the direct of the invariant spire model.
  %%
  %% FUNCTION CALL
  %%
  %% data = directInvariant(sky, bands, Hdirect250, Hdirect360, Hdirect520,
  %% index250, index360, index520, Nalpha, Nbeta, Norder)
  %%
  %% where data is a 3 cells of cells that contains the scans of 250,
  %% 360 and 520 respectively. The scans are cells with a 2D tab in each
  %% cell. The tab is bolometer in line and time in column.
  %%
  %% PARAMETERS
  %%
  %% sky : the input sky. This must be a tab of Nalpha x Nbeta x
  %% (3*Norder). The third dim must be arranged in order 0 of 250, then
  %% 360 then 520 first. After there is order 1 of 250, then 360 then
  %% 520. until Norder.
  %%
  %% bands : is vector of dim 3. If bands(1) equal to 1, the compute for
  %% 250. bands(2) and bands(3) for 360 and 520 respectively.
  %%
  %% Hdirect250 (360/520) : a tab of Nalpha, Nbeta, Norder, Nspeed that
  %% contains the transfert function for 250 (360/520). Nalpha and Nbeta
  %% are the number of pixel in alpha and beta, respectively. Norder is
  %% the number of decomposition in lambda. Nspeed is the number of
  %% speed (typicaly four).
  %%
  %% !!!!!!!!!!!!!!!!!!!!!!!
  %%
  %% !!! IMPORTANT !!  The transfert function must be compute by
  %% respecting the zero convention of the fft2 (when you do Hdirect =
  %% fft2(H), zero of H must be in (1,1)). In addition to avoid phase
  %% problem, H must be zero padded and of size Nalpha x Nbeta. A
  %% solution is to used ri2fourier provided.
  %%
  %% !!!!!!!!!!!!!!!!!!!!!!!
  %%
  %% index250 (360/520) : are the index corresponding to the position in
  %% (alpha, beta) when the data as been aquired for 250 (350/520). This
  %% must be index so use potentialy the matlab function sub2ind.
  %%
  %% Nalpha, Nbeta, Norder : the number of pixel in alpha, beta and the
  %% number of order for lambda decomposition
  %%
  %% FUNCTION CALL
  %%
  %% data = directInvariant(sky, bands, Hdirect250, Hdirect360, Hdirect520,
  %% index250, index360, index520, Nalpha, Nbeta, Norder)

  x_tmp = zeros(Nalpha, Nbeta, Nspeed);
  
  %% Make a convolution for each order for each speed
  for iorder = 0:Norder-1

    for ispeed = 1:Nspeed
      
      %% The convolution
      iorder_s = conv2(sky(:,:,iorder+1),Hdirect(:,:, iorder+1, ...
                                                 ispeed), 'same');
      %% Sum of order (but depend on speed)
      x_tmp(:,:,ispeed) = x_tmp(:,:,ispeed) + iorder_s;
    end

  end

  %% init
  scans = cell(1,Nscan);
  
  %% Prelevement and rearanging
  for iscan = 1:Nscan
    
    %% identify the PSF with the speed
    x_tmp_tmp = x_tmp(:,:,uspeed(1,:) == theSpeeds(1,1,iscan));
    %% Extract the observed pixel for that scan from tmp
    observed = x_tmp_tmp(index{iscan});
    %% This order because index_..._... is order in all sample point
    %% for one bolo, the all sample for the second bolo. All in on
    %% line. You must keep this reshape for coef multiplication in
    %% transpose. Or do the reshape in transpose. Or reshape the
    %% coefs... :)
    scans(iscan) = {reshape(observed, numel(index{iscan})/Nbolo, [])};
  end
  
  data = scans;
  
end
