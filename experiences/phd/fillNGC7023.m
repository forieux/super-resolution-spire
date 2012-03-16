clear all

placemount = '/espace/orieux/results/'
expname = 'ngc7023'
system(['mkdir -p ',placemount,expname]);

disp('Init')

load /espace/orieux/realData/data
load /espace/orieux/realData/pointing

pointing250 = pointing{1};
pointing360 = pointing{1};
pointing520 = pointing{1};

path(path,'/home/orieux/matlab/spire/trunk')
path(path,'/home/orieux/matlab/spire/trunk/libspire')
%path(path,'/home/orieux/matlab/spire/trunk/simulator')
path(path,'/home/orieux/matlab/spire/trunk/utils')
path(path,'/home/seismic/matlab/optimgpi')

for iscan = 1:10
    p = pointing250{iscan};
    p(1,:) = rad2arcsec(deg2rad(p(1,:)).*cos(7.5));
    p(2,:) = rad2arcsec(deg2rad(p(2,:)));
    pointing250(iscan) = {p};

    p = pointing360{iscan};
    p(1,:) = rad2arcsec(deg2rad(p(1,:)).*cos(7.5));
    p(2,:) = rad2arcsec(deg2rad(p(2,:)));
    pointing360(iscan) = {p};
    
    p = pointing520{iscan};
    p(1,:) = rad2arcsec(deg2rad(p(1,:)).*cos(7.5));
    p(2,:) = rad2arcsec(deg2rad(p(2,:)));
    pointing520(iscan) = {p};
end

paramsInstrument
paramsObservation
disp('The speed vector is not good !')
paramsSkyNGC7023

% color = ['b' 'g' 'r' 'c' 'm' 'y'];

% clf
% for iscan = 1:10
%     p = pointing250{iscan};
%     plot(p(1,:),p(2,:),['.',color(mod(iscan,5)+1)])
%     hold on
% end
% axis equal

%% Precomputation
[alpha, beta, pointing250, pointing360, pointing520, index250, index360, ...
 index520, coefs250, coefs360, coefs520] = computeIndex(alpha, beta, ...
                                                  pointing250, pointing360, ...
                                                  pointing520, alpha_step, ...
                                                  beta_step);

% clf
% c = coefs250{1};
% for iscan = 2:10
%     c = c + coefs250{iscan};
% end
% imagesc(c)
% axis equal

%%

disp('Dirty map')
[alpha_mm, beta_mm, pointing250_mm, pointing360_mm, pointing520_mm, ...
 index250_mm, index360_mm, index520_mm, coefs250_mm, coefs360_mm, coefs520_mm] ...
    = computeIndex(alpha, beta, pointing250, pointing360, pointing520, 6, 6);

bands = [1 0 0];

coaddCirrus = dirtymap(data, bands, index250_mm, index360_mm, index520_mm, ...
                       coefs250_mm, coefs360_mm, coefs520_mm, ...
                       length(unique(alpha_mm)) , length(unique(beta_mm)));

data250 = data{1};
a = data250{1};
for iscan = 2:10
    b = data250{iscan};
    a = cat(1,a,b);
end
m = mean(a(1,:));

for iscan = 1:10
    a = data250{iscan};
    c = a(1,:);
    c = repmat(c,size(a,1),1);
    a = a./c*m;
    data250(iscan) = {a};
end
data(1) = {data250};

coaddCirrus2 = dirtymap(data, bands, index250_mm, index360_mm, index520_mm, ...
                        coefs250_mm, coefs360_mm, coefs520_mm, ...
                        length(unique(alpha_mm)) , length(unique(beta_mm)));
 
trou = dirtymap(data, bands, index250, index360, index520, coefs250, ...
                coefs360, coefs520, length(unique(alpha)), ...
                length(unique(beta)), '0');

clear a b c p data250
clear index250_mm index360_mm index520_mm coefs250_mm coefs360_mm coefs520_mm

% subplot(1,2,1)
% imagesc(alpha, beta, coaddCirrus(:,:,1));
% axis equal
% subplot(1,2,2)
% imagesc(alpha, beta, coaddCirrus2(:,:,1) + 2*abs(min2(coaddCirrus2(:,:,1))));
% axis equal

%%

disp('Precalculs RI')
precalculsInvariantRI

%% Data simulation
std250 = (1e-3); gammaB250 = 1/std250^2;
std360 = std250; gammaB360 = 1/std360^2;
std520 = std360; gammaB520 = 1/std520^2;

std = [std250 std360 std520];
gammaB = [gammaB250, gammaB360, gammaB520];


%% GPAC Options
OPS(1)	 = 1;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 1e-3;	% Seuil d'arret sur x
OPS(3)	 = 1e-3;	% Seuil d'arret sur f
OPS(4)	 = 10;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 200;		% Nbre maximum d'itérations
OPS(15)	 = 1e-15;	% Pas minimum (MuEps)
OPS(18)	 = 1e-12;		% Pas a l'origine
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
hypers(2,:,1) = [1 1 1]*1e4;

%% Init
%%% The dirtymap filled with mean
init = dirtymap(data, [1 0 0], index250, index360, index520, coefs250, ...
                coefs360, coefs520, Nalpha, Nbeta, 'M');

%%% Gain adaptation
init(:,:,1) = init(:,:,1)./Hrond250(1);
init(:,:,2) = init(:,:,2)./Hrond360(1);
init(:,:,3) = init(:,:,3)./Hrond520(1);

%%% Prior mean
%%%% No prior mean here (only usefull for sampling)
objectMean = zeros(size(init));

bands = [1 0 0];

OPS(18)	 = 2e-12;		% Pas a l'origine

T0 = cputime;

disp('Opt')

Hrond520 = 0;
Hrond360 = 0;
index360 = 0;
index520 = 0;

T0 = cputime;
[fillmap SortieOPS Histo] = gpac('calcDirtyCrit',init,OPS,'calcDirtyGrad', ...
                                 data, hypers, bands, index250, index360, ...
                                 index520, regOp, Nalpha, Nbeta, 'fig', 1000);
DQ = cputime-T0;
fprintf('\n\n\t Durée = %5.0f secondes i.e. %2.1f minutes \n\n',DQ,DQ/60)

name = [placemount,expname,'/optim',num2str(hypers(2,1))]
save(name, 'fillmap', 'OPS', 'SortieOPS', 'Histo', 'hypers', 'init')


%%

[alpha_mm, beta_mm,  pointing250_mm, pointing360_mm, pointing520_mm, ...
 index250_mm, index360_mm, index520_mm, coefs250_mm, coefs360_mm, coefs520_mm] ...
    = computeIndex(alpha, beta, pointing250, pointing360, pointing520, 6, 6);

coaddNGC = dirtymap(data, bands, index250_mm, index360_mm, index520_mm, ...
                    coefs250_mm, coefs360_mm, coefs520_mm, ...
                    length(unique(alpha_mm)) , length(unique(beta_mm)));

name = [placemount,expname,'/coadd']
save(name, 'coaddNGC', 'alpha_mm', 'beta_mm')

%%

ligne = 820;
ligne_mm = min(find(alpha_mm > alpha(820)));

figure(2)
clf
plot(beta, fillmap(ligne,:,1),'r');
hold on
plot(beta_mm, coaddNGC(ligne_mm,:,1))

                
