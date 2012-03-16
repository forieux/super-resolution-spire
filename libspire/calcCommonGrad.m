function grad = calcCommonGrad(commonSig, output, data, dataBlind, ...
                               hypers, bands, indexCommon, ...
                               psdCommon, commonSigMean, varargin)
  %% CALCCOMMONGRAD - Compute the gradient of the quadratic criterion at
  %% object
  %%
  %% grad = calcCommonGrad(commonSig, output, data, dataBlind, hypers,
  %% bands, indexCommon, psdCommon, commonSigMean)
  %%
  %% calcCommonGrad(..., 'fig', NUMFIG) plot in addition some element on
  %% figure NUMFIG
  %%
  %% INPUT PARAMETERS
  %%
  %% commonSig -- current value of common signal
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
  %% is the index for the scan. These index are the element of commonSig
  %% to take for data adequation. They must correspond to the data that
  %% are used to estimate the sky. So data(index(1)) = H(index(1))x +
  %% commonSig(indexCommon(1))...
  %%
  %% psdCommon -- the prior psd for the commonSig. Equivalent of regOps
  %% but for commonSig.
  %%
  %% commonSigMean -- the mean of the commonSig law.  Must be an
  %% TotalTimeSample x 3.
  %%
  %% FUNCTION CALL
  %%
  %% grad = calcCommonGrad(commonSig, output, data, dataBlind, hypers,
  %% bands, indexCommon, psdCommon, commonSigMean)
    
    %% Extras options
    Nbaseoption = 9;
    if nargin > Nbaseoption
        
        for iargin = 1:2:nargin-Nbaseoption
            if strcmp(varargin{iargin},'fig')
                numfig = varargin{iargin+1};
            else
                disp('Unknown option')
            end
        end
        
    end    
    
    paramsInstrument
    paramsObservation
    
    compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);
    
    output250 = output{1}; output360 = output{2}; output520 = output{3};
    
    data250 = data{1}; data360 = data{2}; data520 = data{3};
    
    %% Compute the diff between model reproduction and data
    error250 = cell(size(data250));
    error360 = cell(size(data360));
    error520 = cell(size(data520));
    
    %% 'e = y - Hx'
    if compute250
        error250 = calcDataError(output250, data250, commonSig(:,1), ...
                                 indexCommon{1}, Nbolo250);
    end
    if compute360
        error360 = calcDataError(output360, data360, commonSig(:,2), ...
                                 indexCommon{2}, Nbolo360);
    end
    if compute520
        error520 = calcDataError(output520, data520, commonSig(:,3), ...
                                 indexCommon{3}, Nbolo520);
    end

    %% Reformat the diff 'e' for transpose function
    
    %% And now for common signal. Very easy
    commonAdeq = zeros(size(commonSig));
    if compute250
        commonAdeq(:,1) = 2*hypers(1,1)*transposeCommon(error250, ...
                                                        commonSig(:,1), ...
                                                        indexCommon{1}, ...
                                                        {dataBlind{1:2}});
    end
    if compute360
        commonAdeq(:,2) = 2*hypers(1,2)*transposeCommon(error360, ...
                                                        commonSig(:,2), ...
                                                        indexCommon{2}, ...
                                                        {dataBlind{3:4}});
    end
    if compute520
        commonAdeq(:,3) = 2*hypers(1,3)*transposeCommon(error520, ...
                                                        commonSig(:,3), ...
                                                        indexCommon{3}, ...
                                                        {dataBlind{5:6}});
    end
    
    %% Regularization
    regGrad = zeros(size(commonAdeq)); 
    if compute250
        cr = myfft(commonSig(:,1) - commonSigMean(:,1));
        regGrad(:,1) =2*hypers(2,1)*real(myifft(cr./psdCommon));
    end

    if compute360
        cr = myfft(commonSig(:,2) - commonSigMean(:,2));
        regGrad(:,2) = 2*hypers(2,2)*real(myifft(cr./psdCommon));
    end
    
    if compute520
        cr = myfft(commonSig(:,3) - commonSigMean(:,3));
        regGrad(:,3) = 2*hypers(2,3)*real(myifft(cr./psdCommon));
    end
    
    %% Full gradient
    grad = commonAdeq + regGrad;
    %% grad = commonAdeq;
    
    %% Plotting
    if exist('numfig','var')

        sfigure(numfig);

        clf
        subplot(1,2,1)
        plot(commonSig(:,2));
        title('360(0)')
        
        subplot(1,2,2)
        plot(grad(:,2));
        title('\Delta 360(0)', 'interpreter', 'tex')

        drawnow
        
    end

end

function error = calcDataError(output, data, commonSig, indexCommon, Nbolo)
%% CALCDATAERROR - Compute the error between data and reproduction
%% 
%% error = calcDataError(output, data, commonSig, indexCommon, Nbolo)
%% 
%% compute the error of reproduction for one array. error is in the format of
%% transpose invariant intput. See it's help (it is a cell of scan).
%% 
%% INPUT PARAMETERS
%% 
%% output -- the data reproduction for these array. Must be a cell of scan
%% with in each cell the data of all the bolometer. See directInvariant help
%% 
%% data -- the data of these array in the same format
%% 
%% commonSig -- the common signal for these array. A vector of N samples
%% where N is the total number of sample.
%% 
%% indexCommon -- the matrix T. A tab of Nscan x Msample where M is the
%% number of samples of data (and also output) for one scan.
%% 
%% Nbolo -- the number of bolometer for these scan
    
  paramsInstrument
  paramsObservation
  
  error = cell(1,N_scan_total);
  
  for iscan = 1:N_scan_total
    
    %% Duplication of the common signal of these scan for each
    %% bolometer. indexCommon is T so compute Tc for these scan and
    %% repmat for reproduction for each bolo.
    common_scan = transpose(repmat(commonSig(indexCommon(iscan,:)), ...
                     [], Nbolo));
    
    common_scan = reshape(common_scan, [], Nbolo);
    
    %% (Hx + Tc - y)
    error(iscan) = {output{iscan} + common_scan - data{iscan}};
    
  end
  
end

function commonAdeq = transposeCommon(error, commonSig, indexCommon, ...
                                      dataBlind)
  %% For one array
  
  paramsInstrument
  paramsObservation
  
  commonAdeq = zeros(size(dataBlind{1}));
  
  for iscan = 1:N_scan_total
    %% Sum of all the error on bolometer dim with sum(error,1). Add this
    %% error to previous scan
    commonAdeq(indexCommon(iscan,:)) = ...
        commonAdeq(indexCommon(iscan,:)) + sum(error{iscan},2)';
  end
  
  commonAdeq = commonAdeq + (commonSig' - dataBlind{1});
  commonAdeq = commonAdeq + (commonSig' - dataBlind{2});
  
end

