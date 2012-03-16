function data = directDirty(sky, index, Nalpha, Nbeta, Nbolo, Nscan)
  %% DIRECTDIRTY - Direct computation of the pixel prelevment
  %%
  %% data = directInvariant(sky, bands, index250, index360, index520,
  %% Nalpha, Nbeta)
  %%
  %% where Data is a 3 cells of cells that contains the scans of 250,
  %% 360 and 520 respectively
  %%
  %% PARAMETERS
  %%
  %% sky -- the input sky. This must be a tab of Nalpha x Nbeta x 3
  %% arranged in 250, then 360 then 520.
  %%
  %%
  %% index250 (360/520) -- are the index corresponding to the position
  %% in (alpha, beta) when the data as been aquired for 250 (350/520).
  %% This must be index so use potentialy the matlab function sub2ind.
  %%
  %% Nalpha, Nbeta -- the number of pixel in alpha, beta.
  %%
  %% FUNCTION CALL
  %%
  %% data = directInvariant(sky, bands, index250, index360, index520,
  %% Nalpha, Nbeta)

    data = cell(1,Nscan);
    
    %Prelevement and rearanging
    for iscan = 1:Nscan
      observed = sky(index{iscan});
      %% This order because index_..._... is order in all sample point for one
      %% bolo, the all sample for the second bolo. All in on line. You
      %% must keep this reshape for coef multiplication in transpose. Or
      %% do the reshape in transpose. Or reshape the coefs... :)
      data(iscan) = {reshape(observed, numel(index{iscan})/Nbolo, [])};
    end

end
