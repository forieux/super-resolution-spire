function [xEap gbChain gChain] = wiener0NS(data, bands, Hrond250, Hrond360, ...
                                           Hrond520, regOp, Nalpha, Nbeta, ...
                                           alphaBound, betaBound, criterion, ...
                                           burnin, maxIter, varargin)
%% WIENER0NS - Wiener/Hunt filter for order zero with hyperparameter
%% estimation
%%
%% This algorithm is a matlab implementation of the method describes in the
%% paper " Unsupervised and myopic deconvolution with a Bayesian approach "
%% by F. Orieux, J.-F. Giovannelli and T. Rodet
%% 
%% FUNCTION CALL
%% 
%% [xEap gbChain gChain] = wiener0NS(data, bands, Hrond250, Hrond360,
%% Hrond520, regOp, Nalpha, Nbeta, alphaBound, betaBound, criterion, burnin,
%% maxIter)
%% 
%% wiener0NS(..., 'fig', figNumber) plot interesting things
%% 
%% wiener0NS(..., 'var', 0 or 1) compute the empirical variance for each
%% pixel if 1 (default 0). Not activated for the moment.
%% 
%% wiener0NS(..., 'bound', 0 or 1) compute the regularity mesure inside the
%% bound provided in paramter if 1. Default is 0.
%% 
%% OUTPUT
%% 
%% xEap -- the posterior mean of the image
%% 
%% gbChain -- the 3xNsample chain of the noise precision
%% 
%% gChain -- the Ncorrelation x 3 x Nsample of the image precision
%%
%% PARAMETERS
%% 
%% data -- the NORMALIZED fourrier transform of retroprojected of data on sky
%% space, with a mean by the number of time the pixel as been observed for
%% each. This must be a 2D TAB of Nalpha x Nbeta. Here Norder is supposed to
%% be 1 (only order 0).
%%    
%% Hrond250 (360/520) -- a tab of Nalpha, Nbeta, Nspeed that contains the
%% DIRECT (the conjugate is automaticlly use) transfert function for 250
%% (360/520). Nalpha and Nbeta are the number of pixel in alpha and beta,
%% respectively. Nspeed is the number of speed (typicaly four).
%%
%% regOp -- a Nalpha x Nbeta x N tab that contains the N regularization
%% operators (a diff operator, a mean operator etc...)
%%
%% Nalpha, Nbeta -- the number of pixel in alpha ans beta.
%%
%% alphaBound, betaBound -- the range in matlab index over which the mesure
%% of regularity must be done. if 1 and end it is over all the image.
%% 
%% criterion -- if the difference between two successive mean is less than
%% this value, stop the algorithm.
%% 
%% burnin -- number of iteration to remove at the beginning of the chain to
%% compute the mean of the image.
%% 
%% maxIter -- maximum number of iteration.
%%
%% FUNCTION CALL
%% 
%% [xEap gbChain gChain] = wiener0NS(data, bands, Hrond250, Hrond360,
%% Hrond520, regOp, Nalpha, Nbeta, alphaBound, betaBound, criterion, burnin,
%% maxIter, varargin)

  fprintf('Initialisation..\n')
  
  bound = 0;
  varianceSwitch = 0;
  
  %% Extras options
    Nbaseoption = 13;
    if nargin > Nbaseoption
        
        for iargin = 1:2:nargin-Nbaseoption
            if varargin{iargin} == 'fig'
                numfig = varargin{iargin+1};
            elseif varargin{iargin} == 'bnd'
                bound = varargin{iargin+1};
            elseif varargin{iargin} == 'var'
                varianceSwitch
            else
                disp('Unknown option')
            end
        end
        
    end    
    
    compute250 = bands(1);
    compute360 = bands(2);
    compute520 = bands(3);
    
    %% For time computation estimation
    tLoop = 0;
    meanTime = 0;
    
    %% The mean of the object
    circXeap = zeros(Nalpha, Nbeta, 3);
    %% The previous computed mean in the iterative loop
    previousCircXeap = zeros(Nalpha, Nbeta, 3);
    
    %% The mean of x^2 to compute the variance
    circX2eap = zeros(Nalpha, Nbeta, 3);
    %% The previous computed mean of x^2
    previousX2eap = zeros(Nalpha, Nbeta, 3);
    
    %% Difference between two succesive mean
    delta = NaN;
    
    %% Initial state of the chain
    circXsample = zeros(Nalpha, Nbeta, 3);

    gbChain = zeros(3,burnin);
    gbChain(:,1) = [1 1 1]*(1e-3)^2;
    
    gChain = zeros(size(regOp,3),3,burnin);
    gChain(:,:,1) = 1e-0*ones(size(regOp,3),3);
    
    %% Parameter of the hyperparameter law. The following value correspond
    %% to Jeffery's prior for the hyper parameter.
    alphaB = 0;
    betaBbar = 0; % = 1/betaB
    
    alpha0 = 0;
    beta0bar = 0; % = 1/beta0
    
    %% Gibbs sampling
    fprintf('Start of sampling...\n')
    
    totalTime = cputime;
    
    for iteration = 2:maxIter
        
        %% For iteration time
        t0 = cputime;
        
        
        precisions = [gbChain(:,iteration - 1) gChain(:,:,iteration - 1)']';
        inverseVariance = hessian(precisions, bands, Hrond250, Hrond360, ...
                                  Hrond520, regOp, Nalpha, Nbeta, 0, 0)./2;
        
        if compute250

            if bound
                [Xsample gbSample gSample] = sampleBandBounded(data(:,:,1), ...
                                                               inverseVariance(:,:,1), ...
                                                               gbChain(1,iteration-1), ...
                                                               Hrond250, regOp, ...
                                                               alphaB, betaBbar, ...
                                                               alpha0, beta0bar, ...
                                                               Nalpha, Nbeta, ...
                                                               alphaBound(:,1), ...
                                                               betaBound(:,1));
            else
                [Xsample gbSample gSample] = sampleBand(data(:,:,1), ...
                                                        inverseVariance(:,:,1), ...
                                                        gbChain(1,iteration-1), ...
                                                        Hrond250, regOp, ...
                                                        alphaB, betaBbar, ...
                                                        alpha0, beta0bar, ...
                                                        Nalpha, Nbeta);
            end
            
            circXsample(:,:,1) = Xsample;
            gbChain(1,iteration) = gbSample;
            gChain(:,1,iteration) = gSample;
            
       end
        
       if compute360
           
           if bound
               [Xsample gbSample gSample] = sampleBandBounded(data(:,:,2), ...
                                                              inverseVariance(:,:,2), ...
                                                              gbChain(2,iteration-1), ...
                                                              Hrond360, regOp, ...
                                                              alphaB, betaBbar, ...
                                                              alpha0, beta0bar, ...
                                                              Nalpha, Nbeta, ...
                                                              alphaBound(:,2), ...
                                                              betaBound(:,2));
           else
               [Xsample gbSample gSample] = sampleBand(data(:,:,2), ...
                                                       inverseVariance(:,:,2), ...
                                                       gbChain(2,iteration-1), ...
                                                       Hrond360, regOp, ...
                                                       alphaB, betaBbar, ...
                                                       alpha0, beta0bar, ...
                                                       Nalpha, Nbeta);
           end

           gSample
           
           circXsample(:,:,2) = Xsample;
           gbChain(2,iteration) = gbSample;
           gChain(:,2,iteration) = gSample;
           
       end
       
       if compute520

           if bound
               [Xsample gbSample gSample] = sampleBandBounded(data(:,:,3), ...
                                                              inverseVariance(:,:,3), ...
                                                              gbChain(3,iteration-1), ...
                                                              Hrond520, regOp, ...
                                                              alphaB, betaBbar, ...
                                                              alpha0, beta0bar, ...
                                                              Nalpha, Nbeta, ...
                                                              alphaBound(:,3), ...
                                                              betaBound(:,3));
           else
               [Xsample gbSample gSample] = sampleBand(data(:,:,3), ...
                                                       inverseVariance(:,:,3), ...
                                                       gbChain(3,iteration-1), ...
                                                       Hrond520, regOp, ...
                                                       alphaB, betaBbar, ...
                                                       alpha0, beta0bar, ...
                                                       Nalpha, ...
                                                       Nbeta);
           end
           
           circXsample(:,:,3) = Xsample;
           gbChain(3,iteration) = gbSample;
           gChain(:,3,iteration) = gSample;
           
       end 
       
       %% ======================================
       %% Current Empirical Mean eq. 31
       
       if iteration > burnin
           
           circXeap = previousCircXeap + circXsample;
           
           %% if varianceSwitch == 1
           %%     X2eap = previousX2eap + real(uifft2(circXsample)).^2;
           %%     previousX2eap = X2eap;
           %% end
           
           if iteration > (burnin + 1)
               norm = sum(abs(circXeap(:) / (iteration - burnin)));
               current = circXeap / (iteration - burnin);
               previous = previousCircXeap / (iteration - burnin - 1);
               
               delta = sum(abs(current(:) - previous(:))) / norm;
           end
            
           previousCircXeap = circXeap;
           
       end
       
       %% ======================================
       %% Plotting
              
       if exist('numfig')
           
           sfigure(numfig);
           set(numfig,'defaultaxesfontsize',11);
           set(numfig,'defaulttextfontsize',11);
           
           clf
           if compute250
               subplot(2,3,1)
               imagesc(real(uifft2(circXsample(:,:,1))))
               setim; axis off
               title('PSW 250')

               subplot(2,3,4)
               plot(squeeze(gbChain(1,:)))
               title('\gamma_b','interpreter','tex')
               
               subplot(2,3,5)
               plot(squeeze(gChain(1,1,:)))
               title('\gamma_\alpha','interpreter','tex')
               
           end
           
           if compute360
               subplot(2,3,2)
               imagesc(real(uifft2(circXsample(:,:,2))))
               setim; axis off
               title('PMW 360')

               subplot(2,3,4); hold on
               plot(squeeze(gbChain(2,:)),'k')
               title('\gamma_b','interpreter','tex')
               
               subplot(2,3,5); hold on
               plot(squeeze(gChain(1,2,:)),'k')
               title('\gamma_\alpha','interpreter','tex')
                              
           end
           
           if compute520
               subplot(2,3,3)
               imagesc(real(uifft2(circXsample(:,:,3))))
               setim; axis off
               title('PLW 520')

               subplot(2,3,4); hold on
               plot(squeeze(gbChain(3,:)),'r')
               title('\gamma_b','interpreter','tex')

               subplot(2,3,5); hold on
               plot(squeeze(gChain(1,3,:)),'r')
               title('\gamma_\alpha','interpreter','tex')

           end 
%            subplot(2,3,4);
%            legend('PSW','PMW','PLW')
%            subplot(2,3,5); 
%            legend('PSW','PMW','PLW')
           
           drawnow
       
       end       
       
       %% ======================================
       %% Stop of the algorithme
       
        if delta < criterion
            
            fprintf(['\n\n\t !! End of sampling by criterion !! Enlapsed time : %.1f sec\n\n'], cputime - totalTime)
            
            circXeap = circXeap / (iteration - burnin) ;
            xEap = real(uifft2(circXeap));
            
            %% if varianceSwitch == 1
            %%     X2eap = X2eap / (iteration - burnin) ;
            %%     xVap = X2eap - xEap.^2;
            %% end
            
            return
            
        end
        
        %% ======================================
        %% loop time computation
        
        tLoop = tLoop + cputime - t0;
        meanTime = tLoop / iteration;
        
        fprintf(['\tIteration %i / %i // e = %.3e\n'],iteration, maxIter, ...
                delta)

        fprintf(['\tLoop time = %.3f / Time left = %.1f sec or %.1f mins\n'], ...
                meanTime, (maxIter - iteration) * meanTime, ((maxIter - ...
                                                          iteration) * ...
                                                          meanTime)/60)     
    end
    
    fprintf(['\nEnd of sampling by maximum iteration !! Enlapsed time : %.1f sec\n'], cputime - totalTime)
    
    %% Empirical mean = EAP
    circXeap = circXeap / (maxIter - burnin) ;
    xEap = real(uifft2(circXeap));
    
    %% if varianceSwitch == 1
    %%     X2eap = X2eap / (maxIter - burnin) ;
    %%     xVap = X2eap - xEap.^2;
    %% end
    
end


function [circXsample gbSample gSample] = sampleBandBounded(data, ...
                                                      inverseVariance, gb, ...
                                                      Hrond, regOp, alphaB, ...
                                                      betaBbar, alphaX, ...
                                                      betaXbar, Nalpha, ...
                                                      Nbeta, alphaBound, ...
                                                      betaBound)
%% SAMPLEBAND - Sample object and parameter for one array
%% 
%% [circXsample gbSample gSample] = sampleBandBounded(data, inverseVariance,
%% gb, Hrond, regOp, alphaB, betaBbar, alphaX, betaXbar, Nalpha, Nbeta,
%% alphaBound, betaBound)
%% 
%% OUTPUT
%% 
%% circXsample -- a sample of the image in Fouier space
%% 
%% gbSample -- a sample of the noise precision
%% 
%% gSample -- a tab of image precision for each correlation
%% 
%% PARAMETERS
%% 
%% data -- the NORMALIZED fourrier transform of retroprojected of data on sky
%% space, with a mean by the number of time the pixel as been observed for
%% each. This must be a 2D TAB of Nalpha x Nbeta. Here Norder is supposed to
%% be 1 (only order 0).
%% 
%% gb -- the precision for the noise
%%    
%% Hrond -- a tab of Nalpha, Nbeta, Nspeed that contains the DIRECT (the
%% conjugate is automaticlly use) transfert function for one of the
%% array. Nalpha and Nbeta are the number of pixel in alpha and beta,
%% respectively. Nspeed is the number of speed (typicaly four).  
%% 
%% regOp -- a Nalpha x Nbeta x N tab that contains the N regularization
%% operators (a diff operator, a mean operator etc...)
%%
%% alphaB, betaBbar, alphaX, betaXbar -- prior parameter law for the gamma
%% law. B is for noise. X is for the image. bar because of matlab convention
%% betaBar is 1/beta.
%% 
%% Nalpha, Nbeta -- the number of pixel in alpha ans beta.
%% 
%% alphaBound, betaBound -- the range in matlab index over which the mesure
%% of regularity must be done. if 1 and end it is over all the image.
    
    Nalpha2 = alphaBound(2) - alphaBound(1);
    Nbeta2 = betaBound(2) - betaBound(1);

    %% Sampling noise precision
    
    dataAdeq = real(uifft2(data - circXsample.*sum(Hrond,4)./4));
    dataAdeq = dataAdeq(alphaBound(1):alphaBound(2),betaBound(1): ...
                        betaBound(2));
    
    likelihood = sum2(dataAdeq(:).^2);
    gbSample = gamrnd(alphaB + Nalpha2*Nbeta2/2, 1 / (betaBbar + likelihood/2));

    %% ======================================
    
    %% Sampling image precision
    
    for iregOp = 1:size(regOp,3)
        
        priorAdeq = real(conj(circXsample).*regOp(:,:,iregOp).*circXsample);
        
        priorAdeq = real(uifft2(priorAdeq));
        priorAdeq = priorAdeq(alphaBound(1):alphaBound(2),betaBound(1): ...
                              betaBound(2));
        
        energy = sum(abs(priorAdeq(:)));
        
        gSample(iregOp) = gamrnd(alphaX + (Nalpha2*Nbeta2 - 1)/2, 1 / ...
                                 (betaXbar + energy/2));
    end
    
    %% ======================================

    %% Sample of f(circX^k | gb^k-1, g1^k-1, phi^k, y)
        
    %When these parameters are fixed, the law of circX is Gaussian. So to
    %optained a sample, we simply need to compute a correlated white
    %gaussian noise and then add the mean (wich is the regularized least
    %square solution, or wiener solution)
    
    %% Normalized white complex gaussian noise with std of 1
    circXcentered = (randn(Nalpha, Nbeta)*sqrt(1/2) + i*randn(Nalpha, ...
                                                      Nbeta)*sqrt(1/2))./ ...
        sqrt(Nalpha*Nbeta);
    
    circXcentered = circXcentered./sqrt(inverseVariance);
    circXmean = gb*sum(conj(Hrond),4).*data./inverseVariance;
    circXsample = circXcentered + circXmean;
    
    %% ======================================
        
end

function [circXsample gbSample gSample] = sampleBand(data, inverseVariance, ...
                                                     gb, Hrond, regOp, ...
                                                     alphaB, betaBbar, ...
                                                     alphaX, betaXbar, ...
                                                     Nalpha, Nbeta)
%% SAMPLEBAND - Sample object and parameter for current band
%% 
%% Since there is no bound, there is no need to go to and return from Fourier
%% space.
%%
%% [circXsample gbSample gSample] = sampleBandBounded(data, inverseVariance,
%% gb, Hrond, regOp, alphaB, betaBbar, alphaX, betaXbar, Nalpha, Nbeta)
%% 
%% OUTPUT
%% 
%% circXsample -- a sample of the image in Fouier space
%% 
%% gbSample -- a sample of the noise precision
%% 
%% gSample -- a tab of image precision for each correlation
%%
%% PARAMETERS
%% 
%% data -- the NORMALIZED fourrier transform of retroprojected of data on sky
%% space, with a mean by the number of time the pixel as been observed for
%% each. This must be a 2D TAB of Nalpha x Nbeta. Here Norder is supposed to
%% be 1 (only order 0).
%% 
%% gb -- the precision for the noise
%%    
%% Hrond -- a tab of Nalpha, Nbeta, Nspeed that contains the DIRECT (the
%% conjugate is automaticlly use) transfert function for one of the
%% array. Nalpha and Nbeta are the number of pixel in alpha and beta,
%% respectively. Nspeed is the number of speed (typicaly four).  
%% 
%% regOp -- a Nalpha x Nbeta x N tab that contains the N regularization
%% operators (a diff operator, a mean operator etc...)
%%
%% alphaB, betaBbar, alphaX, betaXbar -- prior parameter law for the gamma
%% law. B is for noise. X is for the image. bar because of matlab convention
%% betaBar is 1/beta.
%%
%% Nalpha, Nbeta -- the number of pixel in alpha ans beta.   

    %% Sample of f(circX^k | gb^k-1, g1^k-1, phi^k, y)
        
    %% When these parameters are fixed, the law of circX is Gaussian. So to
    %% optained a sample, we simply need to compute a correlated white
    %% gaussian noise and then add the mean (wich is the regularized least
    %% square solution, or wiener solution)
    
    %% Normalized white complex gaussian noise with std of 1
    circXcentered = (randn(Nalpha, Nbeta)*sqrt(1/2) + i*randn(Nalpha, ...
                                                      Nbeta)*sqrt(1/2));
    
    %% I don't really undestand why but it seems better with :
    circXcentered = (randn(Nalpha, Nbeta)*sqrt(1/2) + i*randn(Nalpha, ...
                                                      Nbeta)*sqrt(1/2))./ ...
        sqrt(Nalpha*Nbeta);

    circXcentered = circXcentered./sqrt(inverseVariance);
    circXmean = gb*sum(conj(Hrond),4).*data./(inverseVariance);
    circXsample = circXcentered + circXmean;
    
    %% ======================================
    
    %% Sampling noise precision
    
    dataAdeq = abs(data - circXsample.*sum(Hrond,4)./4).^2;
    dataAdeq = sum(dataAdeq(:));
    gbSample = gamrnd(alphaB + Nalpha*Nbeta/2, 1 / (betaBbar + dataAdeq/2));

    %% ======================================
    
    %% Sampling image precision
    
    for iregOp = 1:size(regOp,3)
        
        priorAdeq = real(conj(circXsample).*regOp(:,:,iregOp).*circXsample);
        priorAdeq = sum(priorAdeq(:));
        gSample(iregOp) = gamrnd(alphaX + (Nalpha*Nbeta - 1)/2, 1 / (betaXbar ...
                                                          + priorAdeq/2));
    
    end
    
end
