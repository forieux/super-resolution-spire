
expname = 'simpleWiener0';
system(['mkdir -p ',place,expname]);

bands = [1 1 1];

hypers = zeros(numReg+1,3,Norder);

% Noise
hypers(1,:,1) = [gammaB250 gammaB360 gammaB520];
hypers(1,:,2) = [gammaB250 gammaB360 gammaB520];
% Diff for O1 only
hypers(2,:,1) = 1e-4;
hypers(3,:,1) = hypers(2,:,1);

%% Init
xchap = wiener0(cielRondInterp, hypers, bands, Hrond250, Hrond360, Hrond520, ...
                regularizationOperators, index250, index360, index520, ...
                Nalpha, Nbeta);

numfig = 6;
numfig = numfig+1;
figure(numfig); clf
imagesc(alpha, beta, xchap(:,:,2))
axis image; axis xy
colorbar; colormap gray
xlabel('\beta', 'interpreter','tex')
ylabel('\alpha', 'interpreter','tex')
title('Wiener0 360')

numfig = numfig+1;
figure(numfig); clf
plot(beta,xchap(200,:,2))
hold on
plot(beta,mapnan(200,:,2)./Hrond360(1,1,1),'k')
plot(beta,sky(200,:,2),'r')
xlabel('\beta', 'interpreter','tex')
title('Slice wiener0 360')

% Time-stamp: < 03/08(aoû)/2009 21:39 by orieux (simpleWiener.m) >
