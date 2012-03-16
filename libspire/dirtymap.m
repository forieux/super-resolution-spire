function map = dirtymap(scans, index, offsets, Nalpha, Nbeta, Nscan, ...
                        Nbolo)
  %% TRANSPOSEMAP - Projection of data on sky space only with pointing
  %%
  %% Compute the transpose of the data on sky space before instrument
  %% effect by coaddition of data
  %%
  %% FUNCTION CALL
  %%
  %% map = dirtymap(data, index, offsets, Nalpha, Nbeta, Nscan, Nbolo)
  %%
  %% map is a Nalpha x Nbeta
  %%
  %% PARAMETERS
  %%
  %% data : the data. This must be a cell of dim 3 that contains 3 cell,
  %% one for 250, 360 and 520 in that order. Each element of cell of
  %% each band must contains the actual data for one scan. For one scan
  %% the data must be a tab arranged to correspond to the correspond
  %% index (the arragement is not important here, this is your job to be
  %% coherent).
  %%
  %% bands : is vector of dim 3. If Params(1) equal to 1, the transpose
  %% for 250 is computed. Params(2) and Params(3) for 360 and 520
  %% respectively.
  %%
  %% index250 (360/520) : are the index corresponding to the position in
  %% (alpha, beta) when the data as been aquired for 250 (350/520). This
  %% must be index so use potentialy the matlab function sub2ind.
  %%
  %% coefs250 (360/520) : are the number of time a pixel as been
  %% observed during one scan. This are cells of 2D tab (representing an
  %% image), one tab for each scan.
  %%
  %% Nalpha, Nbeta : the number of pixel in alpha and beta.
  %%
  %% object is a tab of Nalpha x Nbeta x 3.

  coefs = zeros(Nalpha, Nbeta);
  map = zeros(Nalpha, Nbeta);
  
  for iscan = 1:Nscan
    
    %% Coadition of the data that point the same thing on the sky. To
    %% avoid looping on bolometer, a 2D of (Nalpha x Nbeta) x Nbolo is
    %% used. The data from each bolo will be on a differente column and a
    %% final sum is done. To avoid excess memory consumption, a sparse
    %% matrix is used.
    
    %% There is Nbolo bolometer. "length(index)/Nbolo" is normally the
    %% number of time sample in that scan. Reshape in one row
    indexBolo = reshape(repmat(1:Nbolo, [numel(index{iscan})/Nbolo 1]), [numel(index{iscan}) 1]);
    
    %% Projection on convolued coefficient space. This operation is
    %% critical. It is supposed that the index are arranged in the
    %% bolometer order, and so time order for one bolometer. This must be
    %% arranged like this for index and scans.
    im = sparse(index{iscan}, indexBolo, scans{iscan} - repmat(offsets, size(scans{iscan},1),1), Nalpha*Nbeta, Nbolo);
    imCoefs = sparse(index{iscan}, indexBolo, 1, Nalpha*Nbeta, Nbolo);
    
    %% mean the bolo data (on dim 2), and reshape like an image.
    map = map + reshape(sum(im,2), Nalpha, Nbeta);
    coefs = coefs + reshape(sum(imCoefs,2), Nalpha, Nbeta);
    
  end
  
  %% Mean of all scan
  map(coefs ~= 0) = map(coefs ~= 0)./coefs(coefs ~=0);
  
end

