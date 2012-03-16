clear all

% where result will be saved
placemount = '/espace/orieux/results/'
expname = 'ngc7023-offsets'
system(['mkdir -p ',placemount,expname]);

disp('Init')

% load data and pointing
load /espace/orieux/realData/data
load /espace/orieux/realData/pointing

pointing250 = pointing{1};
pointing360 = pointing{2};
pointing520 = pointing{3};

% set path, specially gpac
path(path,'/home/orieux/codes/spire/trunk')
path(path,'/home/orieux/codes/spire/trunk/libspire')
%path(path,'/home/orieux/codes/spire/trunk/simulator')
path(path,'/home/orieux/codes/spire/trunk/utils')
path(path,'/home/seismic/matlab/optimgpi')

%% Apply Alex formula to get a square map (abstraction of spherical
%coordinate ?)
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

% set parameters
paramsInstrument
paramsObservation
disp('The speed vector is not good !')
paramsSkyNGC7023

%% Precomputation of index and redondancy
[alpha, beta, pointing250, pointing360, pointing520, index250, index360, ...
 index520, coefs250, coefs360, coefs520] = computeIndex(alpha, beta, ...
                                                  pointing250, pointing360, ...
                                                  pointing520, alpha_step, ...
                                                  beta_step);

% Computation of coadd map. Need also new index and redondancy since the
% resolution is not the same.
disp('Dirty map')
[alpha_mm, beta_mm, pointing250_mm, pointing360_mm, pointing520_mm, ...
 index250_mm, index360_mm, index520_mm, coefs250_mm, coefs360_mm, coefs520_mm] ...
    = computeIndex(alpha, beta, pointing250, pointing360, pointing520, 12, 12);

% Compute for the three array
bands = [1 1 1];

%%% ICI

% coaddition

offsets250 = zeros(1,Nbolo250);
offsets360 = zeros(1,Nbolo360);
offsets520 = zeros(1,Nbolo520);

offsets = {offsets250, offsets360, offsets520};

Nalpha_mm = length(unique(alpha_mm));
Nbeta_mm = length(unique(beta_mm));

disp('start')

for iter = 1:100

    coaddCirrus = dirtymapOffsets(data, bands, index250_mm, index360_mm, ...
                                  index520_mm, coefs250_mm, coefs360_mm, ...
                                  coefs520_mm, offsets, Nalpha_mm, Nbeta_mm);

    dataReproduction = directDirty(coaddCirrus, bands, index250_mm, ...
                                   index360_mm, index520_mm, Nalpha_mm, ...
                                   Nbeta_mm);

    offsets = estimOffsets(data, bands, dataReproduction);

end

imagesc(coaddCirrus(:,:,2))
axis square
print -dpng offstetsPMW

imagesc(coaddCirrus(:,:,1))
axis square
print -dpng offstetsPSW

imagesc(coaddCirrus(:,:,3))
axis square
print -dpng offstetsPLW

%%

