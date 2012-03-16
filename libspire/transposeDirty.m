function map = transposeDirty(data, bands, index250, index360, index520, ...
                              Nalpha, Nbeta)
%% TRANSPOSEMAP - Retroprojection of data on sky space
%%
%% FUNCTION CALL
%% 
%% map = transposeDirty(data, bands, index250, index360, index520, Nalpha,
%% Nbeta)
%% 
%% map is a Nalpha x Nbeta x 3 tab.
%% 
%% PARAMETERS
%% 
%% data -- the data. This must be a cell of dim 3 that contains 3 cell,
%% one for 250, 360 and 520 in that order. Each element of cell of each band
%% must contains the actual data for one scan. For one scan the data must be
%% a tab arranged to correspond to the correspond index (the arragement is
%% not important here, this is your job to be coherent).
%%
%% bands -- is vector of dim 3. If Params(1) equal to 1, the transpose for 250
%% is computed. Params(2) and Params(3) for 360 and 520 respectively.
%%
%% index250 (360/520) -- are the index corresponding to the position in
%% (alpha, beta) when the data as been aquired for 250 (350/520). This must
%% be index so use potentialy the matlab function sub2ind.
%%
%% Nalpha, Nbeta -- the number of pixel in alpha and beta.
    
  %% Init
    
  paramsInstrument
  paramsObservation
    
  compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);
    
  scans250 = data{1}; scans360 = data{2}; scans520 = data{3};
    
  %% Final coaddition and Coadition for one scan
  map = zeros(Nalpha, Nbeta, 3);
    
  for iscan = 1:N_scan_total
        
    if compute250 == 1
      im = scansum(scans250{iscan}, index250{iscan}, Nbolo250, ...
                   Nalpha, Nbeta);
      %% Addition to the previous scan data that where acquired with the
      %% same speed.
      map(:,:,1) = map(:,:,1) + im;
    end
        
    if compute360 == 1
      im = scansum(scans360{iscan}, index360{iscan}, Nbolo360, ...
                   Nalpha, Nbeta);
      map(:,:,2) = map(:,:,2) + im;
    end
        
    if compute520 == 1
      im = scansum(scans520{iscan}, index520{iscan}, Nbolo520, ...
                   Nalpha, Nbeta);
      map(:,:,3) = map(:,:,3) + im;
    end
        
  end
    
end

function im = scansum(scan, index, Nbolo, Nalpha, Nbeta)
%% SCANMEAN - Compute the mean of data per scan
%%             
%% im = scanmean(scan, index, coefs, Nbolo, Nalpha, Nbeta)
%% 
%% return an Nalpha x Nbeta matrix that contains the mean of data in scan
%% pointed by index.
%% 
%% INPUTS ARGUMENTS
%% 
%% scan -- a scan for all bolometer of one matrix. It's arregment must be like index
%% 
%% index -- the index corresponding to each sample in scan
%% 
%% Nbolo -- the number of bolo in this matrix
%% 
%% Nalpha, Nbeta -- the number of pixel in alpha and betal

  %% Coadition of the data that point the same thing on the sky. To
  %% avoid looping on bolometer, a 2D of (Nalpha x Nbeta) x Nbolo is
  %% used. The data from each bolo will be on a differente column and a
  %% final sum is done. To avoid excess memory consumption, a sparse
  %% matrix is used.

  %% There is Nbolo bolometer. "length(index)/Nbolo" is normally the
  %% number of time sample in that scan. Reshape in one row
  indexBolo = reshape(repmat(1:Nbolo, [numel(index)/Nbolo 1]), ...
                      [numel(index) 1]);

  %% Projection on convolued coefficient space. This operation is
  %% critical. It is supposed that the index are arranged in the
  %% bolometer order, and so time order for one bolometer. This must be
  %% arranged like this for index and scans.
  im = sparse(index, indexBolo, scan, Nalpha*Nbeta, Nbolo);
    
  %% sum the bolo data (on dim 2), and reshape like an image.
  im = reshape(sum(im,2), Nalpha, Nbeta);
    
end

