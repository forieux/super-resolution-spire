%% INVERSEDESCENTE - Inversion by gradient descent
placemount = '/espace/orieux/results/'
expname = 'simpleDescenteHessian_360';
system(['mkdir -p ',placemount,expname]);

save([placemount,expname,'/initialization'])

%system('rlog -r libspire/transposeInvariant.m');
%system('rlog -r libspire/calcQuadCrit.m');
%system('rlog -r libspire/calcQuadGrad.m');

%% GPAC Options
OPS(1)	 = 1;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 5e-5;	% Seuil d'arret sur x
OPS(3)	 = 5e-5;	% Seuil d'arret sur f
OPS(4)	 = 10;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 150;		% Nbre maximum d'itérations
OPS(15)	 = 1e-15;	% Pas minimum (MuEps)
OPS(18)	 = 5e-13;		% Pas a l'origine
%OPS(19) = 10;		% Sauve toutes le OPS(19) itérations...
%OPS(20) = 1;		% ... dans un fichier portant un numéro
%OPS(22) = 1;		% Positivité

%% Hyperparameters value
% hyper -- the hyperparamter value tab of 3*Norder + 1 lines and 3 columns,
% one column for each band. The line are ordered like this : the first line
% is noise precision. The lines 2 to 2+(Norder-1) is for difference in
% alpha. From 3+(Norder-1) to 3+2*(Norder-1) is for difference in beta. From
% 4+2*Norder to 4+3*(Norder-1) lines is for the mean. From 5+3*(Norder-1) to
% end is for the norm.

%%% Only on operator in line an column in same time
regOp = circDalpha + circDbeta;

% Noise
hypers = zeros(2,3,1);
hypers(1,:,1) = [gammaB250 gammaB360 gammaB520];
hypers(1,:,2) = [gammaB250 gammaB360 gammaB520];
% Diff for O1 only
hypers(2,:,1) = 3.9e12;

%% Init
%%% The dirtymap filled with mean
init = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

%%% Or load previous computed filled map (see simpleFill)
% load /espace/orieux/results/simpleFill/optim500000 fillmap;
% init = fillmap;

%%% Gain adaptation
init(:,:,1) = init(:,:,1)./Hrond250(1);
init(:,:,2) = init(:,:,2)./Hrond360(1);
init(:,:,3) = init(:,:,3)./Hrond520(1);

%%% Prior mean
%%%% No prior mean here (only usefull for sampling)
objectMean = zeros(size(init));

%% Computation for 360 only
bands = [0 1 0];

OPS(18)	 = 1e-12;		% Pas a l'origine

T0 = cputime;

[xchap SortieOPS Histo] = gpac('calcQuadCrit', init, OPS,'calcQuadGrad', ...
                               data, hypers, bands, Hrond250, Hrond360, ...
                               Hrond520, index250, index360, index520, ...
                               regOp, objectMean, Nalpha, Nbeta, Norder, ...
                               'fig', 1000);

DQ = cputime-T0;
fprintf('\n\n\t Durée = %5.0f secondes i.e. %2.1f minutes \n\n',DQ,DQ/60)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hessian correction
h00 = hessian(hypers, bands, Hrond250, Hrond360, Hrond520, regOp, Nalpha, ...
              Nbeta, 0, 0);

OPS(5)	 = 0; % Strategie 0 == gradient
OPS(7)	 = 1; % min en ligne
OPS(18)	 = 1; % pas à l'origine

T0 = cputime;
[xchapHess SortieOPSHess HistoHess] = gpac('calcQuadCrit', init, ...
                                           OPS,'calcQuadGrad', data, hypers, ...
                                           bands, Hrond250, Hrond360, ...
                                           Hrond520, index250, index360, ...
                                           index520, regOp, objectMean, ...
                                           Nalpha, Nbeta, Norder, 'fig', ...
                                           1000, 'hes', h00);

DQ = cputime-T0;
fprintf('\n\n\t Durée = %5.0f secondes i.e. %2.1f minutes \n\n',DQ,DQ/60)

name = [placemount,expname,'/compHessian_noHessian_0_360_',num2str(hypers(2,1,1))]
save(name, 'xchapHess', 'SortieOPSHess', 'HistoHess', 'xchap', 'OPS', 'SortieOPS', 'Histo', 'hypers', 'init')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
semilogy(Histo(:,4), Histo(:,1))
hold on
semilogy(HistoHess(:,4), HistoHess(:,1),'r')



