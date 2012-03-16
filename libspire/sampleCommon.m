function [comSample SortieOPS Histo] = sampleCommon(init, options, ...
                                                    output, data, ...
                                                    dataBlind, ...
                                                    hypers, bands, ...
                                                    indexCommon, ...
                                                    psdCommon, varargin)
  %% SAMPLECOMMON - Simulate the conditionnal posterior law of common
  %% signal
  %%
  %% [comSample SortieOPS Histo] = sampleCom(init, options, output,
  %% data, dataBlind, hypers, bands, indexCommon, psdCommon)
  %%
  %% compute a sample of the quadratic conditionnal posterior law of the
  %% common signal by optimization.
  %%
  %% Additionnal output are the results of the optimisation see 'gpac'
  %% help.
  %%
  %% It use the fonction 'calcCommonCrit', 'calcCommonGrad' and 'gpac'
  %% to compute the optimization
  %%
  %% sampleCommon(..., 'fig', NUMFIG) plot in addition some element on
  %% figure NUMFIG
  %%
  %% INPUT PARAMETERS
  %%
  %% commonSig -- current value of common signal
  %%
  %% options -- the option vector for GPAC, see it's help.
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
  %% FUNCTION CALL
  %%
  %% [comSample SortieOPS Histo] = sampleCom(init, options, output,
  %% data, dataBlind, hypers, bands, indexCommon, psdCommon)

    paramsInstrument
    paramsObservation
    
    compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);
    
    gammaB250 = hypers(1,1); gammaB360 = hypers(1,2); gammaB520 = hypers(1,3);
    
    data250 = data{1}; data360 = data{2}; data520 = data{3};
    
    %% Data perturbation
    for iscan = 1:N_scan_total
        if compute250
            data250(iscan) = {data250{iscan} + randn(size(data250{iscan}))/ ...
                              sqrt(gammaB250)};
        end
        if compute360
            data360(iscan) = {data360{iscan} + randn(size(data360{iscan}))/ ...
                              sqrt(gammaB360)};
        end
        if compute520
            data520(iscan) = {data520{iscan} + randn(size(data520{iscan}))/ ...
                              sqrt(gammaB520)};
        end
    end
    
    data = {data250 data360 data520};

    if compute250
        dataBlind(1) = {dataBlind250{1} + randn(size(dataBlind{1}))/ ...
                        sqrt(gammaB250)};
        dataBlind(2) = {dataBlind250{2} + randn(size(dataBlind{2}))/ ...
                        sqrt(gammaB250)};
    end
    if compute360
        dataBlind(3) = {dataBlind250{3} + randn(size(dataBlind{3}))/ ...
                        sqrt(gammaB360)};
        dataBlind(4) = {dataBlind250{4} + randn(size(dataBlind{4}))/ ...
                        sqrt(gammaB360)};
    end
    if compute520
        dataBlind(5) = {dataBlind250{5} + randn(size(dataBlind{5}))/ ...
                        sqrt(gammaB520)};
        dataBlind(6) = {dataBlind250{6} + randn(size(dataBlind{6}))/ ...
                        sqrt(gammaB520)};
    end
        
    
    %% Sky perturbation 
    commonSigMean = zeros(size(dataBlind{1}),3);
    
    if compute250
        commonSigMean(:,1) = randn(size, dataBlind{1},1);
        commonSigMean(:,1) = real(myifft(myfft(commonSigMean(:,1)).* ...
                                         sqrt(psdCommon)))/sqrt(hypers(2,1));
    end
    
    if compute360
        commonSigMean(:,2) = randn(size, dataBlind{1},1);
        commonSigMean(:,2) = real(myifft(myfft(commonSigMean(:,2)).* ...
                                    sqrt(psdCommon)))/sqrt(hypers(2,2));
    end
    
    if compute520
        commonSigMean(:,3) = randn(size, dataBlind{1},1);
        commonSigMean(:,3) = real(myifft(myfft(commonSigMean(:,3)).* ...
                                         sqrt(psdCommon)))/sqrt(hypers(2,3));
    end
    
    %% Optimization
    [comSample SortieOPS Histo] = gpac('calcCommonCrit', init, options, ...
                                       'calcCommonGrad', output, data, ...
                                       dataBlind, hypers, bands, indexCommon, ...
                                       psdCommon, commonSigMean, varargin{:});

end
