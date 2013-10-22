function [skyEap gnChain gxChain offsetsChain] = usmse(init, hypersInit, data, ...
                                          Hrond, index, coefs, ...
                                          offsets, regOps, Nalpha, ...
                                          Nbeta, Norder, Nscan, ...
                                          Nbolo, Nspeed, uspeed, ...
                                          theSpeeds, criterion, ...
                                          burnin, maxIter, ...
                                          optptions, varargin)
  %% USMSE - Unsupervised mean square error
  %%
  %% This algorithm implement an unsupervised (hyper parameters
  %% estimation) minimum regularized mean squarre error. This is done
  %% with a Gibbs algo that sample the image the the hyper parameter in
  %% a recursive way.
  %%
  %% The sampling of the image is done with an optimisation algo.
  %% Suppose the same noise level for all bolometers in an array.
  %%
  %% FUNCTION CALL
  %%
  %%
  %%
  %% PARAMETERS
  %%
  %% init -- the first sky sample.
  %%
  %% hypersInit -- the first hypers-paramters sample. If there are not
  %% zero, ignore the first sky sample init, and use these value as
  %% directly to compute a first sample of the sky. A tab of 2 lines
  %% (noise then sky correlation) and 3 column one for each array.
  %%
  %% data -- are the data in directInvariant output convention. A three
  %% cell one for each band. Each cell contains data for all the scan
  %% one cell for each scan.
  %%
  %% Hrond250 (360/520) -- a tab of Nalpha, Nbeta, Nspeed that contains
  %% the DIRECT (the conjugate is automaticlly use) transfert function
  %% for 250 (360/520). Nalpha and Nbeta are the number of pixel in
  %% alpha and beta, respectively. Nspeed is the number of speed
  %% (typicaly four).
  %%
  %% regOps -- a Nalpha x Nbeta x Norder tab that contains the inverse
  %% variance (a diff operator, a mean operator etc...). The same of the
  %% three array.
  %%
  %% Nalpha, Nbeta, Norder -- the number of pixel in alpha, beta and the
  %% number of wavelength developmenent.
  %%
  %% criterion -- if the difference between two successive mean is less
  %% than this value, stop the algorithm.
  %%
  %% burnin -- number of iteration to remove at the beginning of the
  %% chain to compute the mean of the image.
  %%
  %% maxIter -- maximum number of iteration.
  %%
  %% EXAMPLES
  %%
  %% FUNCTION CALL

  %% Extras options
  Nbaseoption = 20;
  if nargin > Nbaseoption

    for iargin = 1:2:nargin-Nbaseoption
      if strcmp(varargin{iargin},'pla')
        outputDir = varargin{iargin+1};
        [SUCCESS,MESSAGE,MESSAGEID] = mkdir(outputDir,'samples');
        if SUCCESS == 0
          disp(MESSAGE)
        end
      else
        disp('Unknown option')
      end
    end

  end

  %% Number of data
  Ndata = 0;
  for iscan = 1:Nscan
    Ndata = Ndata + numel(data{iscan});
  end

  %% Difference between two succesive mean
  delta = NaN;

  %% Go in Fourier space
  skySample = ufft2(init);
  skyEap = zeros(size(skySample));
  previousSkyEap = zeros(size(skySample));

  %% Precision (inverse variance) of noise and sky

  %% Parameter of the hyperparameter law. The following value correspond
  %% to Jeffery's prior for the hyper parameter. Same prior law for gx
  %% an gn
  alphaB = 0;
  betaBbar = 0; % = 1/betaB

  alphaX = 0;
  betaXbar = 0; % = 1/beta0

  gnSample = hypersInit(1,1);
  gxSample = hypersInit(2,1);

  gnChain = zeros(1,burnin);
  gnChain(1) = gnSample;

  gxChain = zeros(Norder,burnin);
  gxChain(:,1) = gxSample;

  offsetsChain = zeros(length(offsets), burnin);

  %% If the hyper init is not provided
  if gnChain(1,1) == 0
    %% ======================================
    %% Sampling of noise precision

    output = directInvariantF(skySample, Hrond, index, Nalpha, ...
                              Nbeta, Norder, Nbolo, Nspeed, Nscan, ...
                              uspeed, theSpeeds);

    %% ||y - Hx||^2
    adeq = 0;
    for iscan = 1:N_scan_total
      adeq = adeq + sum(sum((data{iscan} - output{iscan}).^2));
    end

    gnSample = gamrnd(alphaB + Ndata/2, 1/(betaBbar + adeq/2));

    gnChain(1) = gnSample;

    %% ======================================
    %% Sampling of the object precision


    for iorder = 1:Norder

      %% gammaI (Fx^dag) Q_I (Fx)
      regularity = ...
          real(conj(skySample(:,:,iorder)).*regOps(:,:,iorder).* ...
               skySample(:,:,iorder));

      %% modified
      gxSample(iorder) = gamrnd(alphaX + (Nalpha*Nbeta - 1)/2, ...
                                1/(betaXbar + sum(regularity(:))/2));

    end
    gxChain(:,1) = gxSample;

  end

  skySample = init;

  offsetsSample = offsets;
  fixedOffset = 1;

  %% Gibbs sampler
  for iteration = 2:maxIter

    %% For iteration time
    t0 = cputime;

    %% ======================================
    %% Sampling of the sky

    %% Simulation by optimization.

    %% hypers -- the hyperparamter value tab of 2 x 3 x Norder. The
    %% first line hypers(1,:,1) is the noise precision (so indep of
    %% order) for each band. The second line is the hypers value for
    %% the object. Column is the band.

    hypers = zeros(2,Norder);
    hypers(1,1) = gnSample;
    hypers(2,:) = gxSample;

    [skySample ooptim] = sampleSky(skySample, optptions, data, hypers, ...
                                   Hrond, index, coefs, offsetsSample, ...
                                   regOps, Nalpha, Nbeta, Norder, ...
                                   Nscan, Nbolo, Nspeed, uspeed, ...
                                   theSpeeds);
    if exist('outputDir','var')
      try
        name = [outputDir,'/samples/sample_',num2str(iteration)];
        save(name, 'skySample', 'gxSample', 'gnSample')
      catch exception
      end
    end

    %% ======================================
    %% Sampling of noise precision

    output = directInvariantF(skySample, Hrond, index, Nalpha, ...
                              Nbeta, Norder, Nbolo, Nspeed, Nscan, ...
                              uspeed, theSpeeds);

    adeq = 0;
    for iscan = 1:Nscan
      adeq = adeq + sum(sum((data{iscan} - output{iscan}).^2));
    end

    gnSample = gamrnd(alphaB + Ndata/2, 1/(betaBbar + adeq/2));
    gnChain(iteration) = gnSample;

    %% ======================================
    %% Sampling of the object precision

    for iorder = 1:Norder

      %% gammaI (Fx^dag) Q_I (Fx)
      regularity = ...
          real(conj(skySample(:,:,iorder)).*regOps(:,:,iorder).* ...
               skySample(:,:,iorder));

      %% modified
      gxSample(iorder) = gamrnd(alphaX + (Nalpha*Nbeta - 1)/2, ...
                                1/(betaXbar + ...
                                   sum(regularity(:))/2));

    end
    gxChain(:,iteration) = gxSample;

    %% ======================================
    %% Sampling of offsets

    output = directInvariantF(skySample, Hrond, index, Nalpha, ...
                              Nbeta, Norder, Nbolo, Nspeed, Nscan, ...
                              uspeed, theSpeeds);

    offsetsMean = estimOffsets(data, output, Nscan);
    offsetsSample = offsetsMean + randn(size(offsetsMean)) ./ sqrt(gnSample);
    offsetsSample(fixedOffset) = offsets(fixedOffset);
    offsetsChain(:,iteration) = offsetsSample;

    %% ======================================
    %% Current Empirical Mean

    if iteration > burnin

      skyEap = previousSkyEap + skySample;

      if iteration > (burnin + 1)
        norm = sum(abs(skyEap(:) / (iteration - burnin)));
        current = skyEap / (iteration - burnin);
        previous = previousSkyEap / (iteration - burnin - 1);

        delta = sum(abs(current(:) - previous(:))) / norm;
      end

      previousSkyEap = skyEap;

    end

    %% ======================================
    %% Plotting

    numfig = 2001;
    if exist('numfig','var')

      sfigure(numfig);
      clf
      subplot(1,2,1)
      imagesc(real(uifft2(skySample(:,:,1))))
      axis image; colormap(hot); axis off

      subplot(1,2,2)
      imagesc(real(uifft2(skyEap(:,:,1))))
      axis image; colormap(hot); axis off

      sfigure(numfig+1);
      clf
      subplot(2,2,1)
      plot(squeeze(gnChain(:)),'.')
      title(['gn = ',num2str(gnSample)])

      subplot(2,2,2)
      plot(squeeze(gxChain(1,:)),'.')
      title(['gx = ',num2str(gxSample)])

      subplot(2,2,3)
      plot(log(squeeze(gnChain(:))),'.')
      title('log(gn)')

      subplot(2,2,4)
      plot(log(squeeze(gxChain(1,:))),'.')
      title('log(gx)')

      drawnow
    end

    %% ======================================
    %% Stop of the algorithme

    if delta < criterion

      skyEap = skyEap / (iteration - burnin) ;
      skyEap = real(uifft2(skyEap));
      return

    end

    tLoop = (cputime - t0)/60;
    leftTime = tLoop*(maxIter - iteration);

    fprintf('Iter %i/%i // e=%.3e // Ltime=%.3f'' / Ttime=%.1f''/%.1f h/%.1f d\n', iteration, maxIter, delta, tLoop, leftTime, leftTime/60, leftTime/(60*24))

  end
  %% ======================================
  %% loop time computation

  %% Empirical mean = EAP
  skyEap = skyEap / (maxIter - burnin);
  skyEap = real(uifft2(skyEap));

end

