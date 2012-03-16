function offsets = estimOffsets(data, dataRepro, Nscan)
  %% ESTIMOFFSETS - Estimation of offset
  %%
  %% The estimation is based on data and data reprodution at current
  %% point and correspond to the mean of the residual.

  offsets = zeros(1,size(data{1},2));
  N = 0;

  for iscan = 1:Nscan
    offsets = offsets + sum(data{iscan} - dataRepro{iscan});
    N = N + size(data{iscan},1);
  end
  offsets = offsets./N;

end
    
