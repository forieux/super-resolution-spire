function object = transposeInvariant2(data, Hdirect, index, Nalpha, ...
                                     Nbeta, Norder, Nbolo, Nspeed, Nscan, ...
                                     uspeed, theSpeeds, varargin)
  %% TRANSPOSEINVARIANT - transpose of the invariant spire model
  %%
  %% This code compute the transpose of the invariant spire model using
  %% the conjugate property
  %%
  %% FUNCTION CALL
  %%
  %% [object] = transposeInvariant(Data, Params, Hdirect250, Hdirect360,
  %% Hdirect520, index250, index360, index520, Nalpha, Nbeta, Norder)
  %%
  %% PARAMETERS
  %%
  %% Data : the data. This must be a cell of dim 3 that contains 3 cell,
  %% one for 250, 360 and 520 in that order. Each element of cell of
  %% each band must contains the actual data for one scan. For one scan
  %% the data must be a tab arranged to correspond to the correspond
  %% index (the arragement is not important here, this is your job to be
  %% coherent).
  %%
  %% Params : is vector of dim 3. If Params(1) equal to 1, the transpose
  %% for 250 is computed. Params(2) and Params(3) for 360 and 520
  %% respectively.
  %%
  %% Hdirect250 (360/520) : a tab of Nalpha, Nbeta, Norder, Nspeed that
  %% contains the DIRECT (the conjugate is automaticlly use) transfert
  %% function for 250 (360/520). Nalpha and Nbeta are the number of
  %% pixel in alpha and beta, respectively. Norder is the number of
  %% decomposition in lambda. Nspeed is the number of speed (typicaly
  %% four).
  %%
  %% !!!!!!!!!!!!!!!!!!!!!!!
  %%
  %% !!! IMPORTANT !!  The transfert function must be compute by
  %% respecting the zero convention of the fft2 (when you do Hdirect =
  %% fft2(H), zero of H must be in (1,1)). In addition to avoid phase
  %% problem, H must be zero padded and of size Nalpha x Nbeta. A
  %% solution is to use the provided ri2fourier.
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
  %% [object] = transposeInvariant(Data, Params, Hdirect250, Hdirect360,
  %% Hdirect520, index250, index360, index520, Nalpha, Nbeta, Norder)

  %% There is two loop in this code, one for the cooadition of the data
  %% on sky (loop on scan) and one for the convolution.

  %% Final coaddition
  x_tmp = zeros(Nalpha, Nbeta, Nspeed);
  
  for iscan = 1:Nscan
    
    %% identify the PSF with the speed
    speed_index = find(uspeed(1,:) == theSpeeds(1,1,iscan));
    
    %% Coadition of the data that point the same thing on the sky. To
    %% avoid looping on bolometer, a 2D of (Nalpha x Nbeta) x Nbolo is
    %% used. The data from each bolo will be on a differente column and
    %% a final sum is done. To avoid excess memory consumption, a sparse
    %% matrix is used.
    
    %% There is 139 bolometer in 250. "length(index250{iscan})/139" is
    %% normally the number of time sample in that scan.
    indexBolo = repmat(1:Nbolo, [length(index{iscan})/Nbolo 1]);
    %% In one row
    indexBolo = indexBolo(:);
    
    %% Projection on convolued coefficient space. This operation is
    %% critical. It is supposed that the index250 are arranged in the
    %% bolometer order, and so time order for one bolometer. This must
    %% be arranged like this for index and scans.
    x_scan = sparse(index{iscan}, indexBolo, data{iscan}, ...
                    Nalpha*Nbeta, Nbolo);
    %% Sum the bolo data (on dim 2), and reshape like an image. Addition
    %% to the previous scan data that where acquired with the same
    %% speed.
    x_tmp(:,:,speed_index) = x_tmp(:,:,speed_index) + ...
        reshape(sum(x_scan,2), Nalpha, Nbeta);

  end

  object = zeros(Nalpha, Nbeta, Norder);

  %% For each speed and each order, to the convolution in Fourier space with
  %% the transpose or the conjugate.
  for ispeed = 1:Nspeed
    
    for iorder = 0:Norder-1
      %% Convolution in Fourier space and addition for one order off all scan even
      %% with different TF.
      
      obj = conv2(x_tmp(:,:,ispeed), ...
                  fliplr(flipud(Hdirect(:,:,iorder+1,ispeed))), 'same');
      object(:,:,iorder*3+1) = object(:,:, iorder*3+1) + obj;
      
    end
  end
end
