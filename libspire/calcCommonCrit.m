function crit = calcCommonCrit(commonSig, output, data, dataBlind, ...
                               hypers, bands, indexCommon, ...
                               psdCommon, commonSigMean, varargin)
  %% CALCCOMMONCRIT - Compute the value of the quadratic criterion value
  %% of object
  %%
  %% val = calcCommonCrit(commonSig, output, data, dataBlind, hypers,
  %% bands, indexCommon, psdCommon, commonSigMean)
  %%
  %% compute the criterion value at object. Data adequation and
  %% regularization are quadratic.
  %%
  %% INPUT PARAMETERS
  %%
  %% commonSig -- current value of common signal. A tab of Nx3 on column
  %% for each array.
  %%
  %% output -- is directInvariant(object, bands, Hrond250, Hrond360,
  %% Hrond520, index250, index360, index520, Nalpha, Nbeta, Norder)
  %%
  %% data -- are the data in directInvariant output convention. A three
  %% cell one for each band. Each cell contains data for all the scan
  %% one cell for each scan.
  %%
  %% dataBlind -- the data of blind bolometer. It is a cell ordered in :
  %% first blind for 250, second for 250, first 360, second 360, first
  %% 520, second 520.
  %%
  %% hypers -- the hypersparamter value tab of (N+1) x 3 x Norder with
  %% the number of line is the number of regularization operator, for
  %% each band in column, and the third dimension for the order. The
  %% first line hypers(1,:,1) is the noise precision (so indep of order)
  %% for each band.
  %%
  %% bands -- a vector of 3 with 1 to compute a band ordered in 250, 360
  %% and 520. Ex : [0 1 0] compute only for 360.
  %%
  %% indexCommon -- (FIXME : adapte to different number of sample per
  %% scan). A cell of 3, one for 250, 360 and 520 respectively. Each
  %% cell contains a Nscan x Nsample tab of index. Line is scan. Column
  %% is the index for the scan. These index are the element of
  %% commonSignal to take for data adequation. They must correspond to
  %% the data that are used to estimate the sky. So data(index(1)) =
  %% H(index(1))x + commonSignal(indexCommon(1))...
  %%
  %% psdCommon -- the prior psd for the commonSig. Equivalent of regOps
  %% but for commonSig.
  %%
  %% commonSigMean -- the mean of the commonSignal law.  Must be an
  %% TotalTimeSample x 3.
  %%
  %% FUNCTION CALL
  %%
  %% val = calcCommonCrit(commonSig, output, data, dataBlind, hypers,
  %% bands, indexCommon, psdCommon, commonSigMean)

  crit = 0;
  
  %% Init
  paramsInstrument
  paramsObservation
  
  compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);    
  
  %% Data adequation
  %% Data reproduction Hx
  
  data250 = data{1}; data360 = data{2}; data520 = data{3};
  output250 = output{1}; output360 = output{2}; output520 = output{3};
  
  %% ||y - Hx||^2
  if compute250
    adeq250 = calcDataAdeq(output250, data250, commonSig(:,1), ...
                           indexCommon{1}, Nbolo250);
    adeq250blind = sum((dataBlind{1} - commonSig(:,1)').^2);
    adeq250blind = adeq250blind + sum((dataBlind{2} - commonSig(:,1)').^2);
    
    crit = hypers(1,1)*(adeq250 + adeq250blind);
  end
  if compute360
    adeq360 = calcDataAdeq(output360, data360, commonSig(:,2), ...
                           indexCommon{2}, Nbolo360);
    adeq360blind = sum((dataBlind{3} - commonSig(:,2)').^2);
    adeq360blind = adeq360blind + sum((dataBlind{4} - commonSig(:,2)').^2);
    
    crit = crit + hypers(1,2)*(adeq360 + adeq360blind);
  end
  if compute520
    adeq520 = calcDataAdeq(output520, data520, commonSig(:,3), ...
                           indexCommon{3}, Nbolo520);
    adeq520blind = sum((dataBlind{5} - commonSig(:,3)').^2);
    adeq520blind = adeq520blind + sum((dataBlind{6} - commonSig(:,3)').^2);
    
    crit = crit + hypers(1,3)*(adeq520 + adeq520blind);
  end
  
  if compute250
    cr = myfft(commonSig(:,1) - commonSigMean(:,1));
    crit = crit + hypers(2,1)*sum(real(myifft(conj(cr).*cr./psdCommon)));
  end
  
  if compute360
    cr = myfft(commonSig(:,2) - commonSigMean(:,2));
    crit = crit + hypers(2,2)*sum(real(myifft(conj(cr).*cr./psdCommon)));
  end
  
  if compute520
    cr = myfft(commonSig(:,3) - commonSigMean(:,3));
    crit = crit + hypers(2,3)*sum(real(myifft(conj(cr).*cr./psdCommon)));
  end

end

function adeq = calcDataAdeq(output, data, commonSig, indexCommon, ...
                             Nbolo)
  %% CALCDATAADEQ - Compute the data adequation for one array
  %%
  %% adeq = calcDataAdeq(output, data, commonSig, indexCommon, Nbolo)
  %%
  %% compute the ||y - Hx - Tc||^2 term of the criterion
  %%
  %% INPUT PARAMETERS
  %%
  %% output -- the data reproduction for these array. Must be a cell of
  %% scan with in each cell the data of all the bolometer. See
  %% directInvariant help
  %%
  %% data -- the data of these array in the same format
  %%
  %% commonSig -- the common signal for these array. A vector of N
  %% samples where N is the total number of sample.
  %%
  %% indexCommon -- the matrix T. A tab of Nscan x Msample where M is
  %% the number of samples of data (and also output) for one scan.
  %%
  %% Nbolo -- the number of bolometer for these scan
    
    adeq = 0;
    for iscan = 1:N_scan_total
        
      %% Duplication of the common signal of these scan for each
      %% bolometer
      common_scan = repmat(commonSig(indexCommon(iscan,:)), [], Nbolo);
      common_scan = reshape(common_scan, [], Nbolo);
        
      %% ||y - Hx - Tc||^2
      adeq = adeq + sum(sum((data{iscan} - output{iscan} - ...
                             common_scan).^2));
   
    end

end

