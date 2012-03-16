function [map offsets] = dirtymapOffsets(scans, offsets, index, Nalpha, Nbeta, Nscan, Nbolo, repeat)
  %% TRANSPOSEMAP - Retroprojection of data on sky space
  %%
  %% This code compute the transpose of the data on sky space before the
  %% duplication of order and the convolution. More exactly, in addition
  %% to the retropprojection, there is a mean instead of a sum of the
  %% pixel that where observed by several data.
  %%
  %% FUNCTION CALL
  %%
  %% map = transposeMap(data, bands, index250, index360, index520,
  %% coefs250, coefs360, coefs520, Nalpha, Nbeta)
  %%
  %% map is a Nalpha x Nbeta x 3 tab.
  %%
  %% or
  %%
  %% map = transposeMap(data, bands, index250, index360, index520,
  %% coefs250, coefs360, coefs520, Nalpha, Nbeta, 'mean')
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

  for iRun = 1:repeat
%     disp(['Estimation of offsets for coaddition (',num2str(iRun),'/',num2str(repeat),')']);
      
    map = dirtymap(scans, index, offsets, Nalpha, Nbeta, Nscan, Nbolo);
    dataRepro = directDirty(map, index, Nalpha, Nbeta, Nbolo, Nscan);
    offsets = estimOffsets(scans, dataRepro, Nscan);
    
  end
  
end
