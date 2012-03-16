clear all

format('long');

load ../data&true/data_fits_image
load ../data&true/pos_fits_image

load mycm

racine = '../SPIRE/';
expname = '13avril/';
system(['mkdir -p ',racine,expname]);
fits_name = 'result_PSW_6arcsec_13avril_pv.fits';
gpac_name = '/gpac_PSW_13avril_fits_pv';

%% set path, specially gpac
path(path,'../')
path(path,'../libspire')
path(path,'../utils')
path(path,'/Users/seismic/matlab/optimgpi')

scan_tot = 10;
pointing250 = pointing{1};
pointing360 = pointing{2};
pointing520 = pointing{3};

temp_ra=[];
temp_dec=[];

%% determining ra_mean and dec_mean
for iscan = 1:scan_tot
  p = pointing250{iscan};
  temp_ra=[temp_ra deg2rad(p(1,:))];
  temp_dec=[temp_dec deg2rad(p(2,:))];
end

raz250=mean(temp_ra);
decz250=mean(temp_dec);

tan_ra=deg2rad(315.47608333333335);
tan_dec=deg2rad(68.15972222222223);

%% Même pixel avant et après projection pour la coordonnée de référence
%% (quelque soit le sytême de coordonnée)
ptot1=[];
ptot2=[];
for iscan = 1:scan_tot
  p = pointing250{iscan};
  ptot1 = union(p(1,:),ptot1);
  ptot2 = union(p(2,:),ptot2);
end
ptot = [ptot1;ptot2];
ptot = deg2rad(ptot);
[xi,eta] = convert_ra_dec_bis(ptot, tan_ra, tan_dec );

%% Projection
for iscan = 1:scan_tot
  p = pointing250{iscan};
  p(2,:) = deg2rad(p(2,:));
  p(1,:) = deg2rad(p(1,:));
  [xi,eta] = convert_ra_dec_bis( p, tan_ra, tan_dec );
  z(1,:) = rad2arcsec(xi);
  z(2,:) = rad2arcsec(eta);
  pointing250(iscan) = {z};
  clear z
end

paramsInstrumentNGC7023
paramsObsNGC7023
paramsSkyNGC7023

bound = 4*sigma_coef*max(band_520) + 4*60; % arcsec
[alpha beta] = computeAxis(pointing250, sky_alpha_period, sky_beta_period, ...
                          bound, N_scan_total);

%% Precomputation of index and redondancy
[alphaCoadd, betaCoadd, pointing250coadd, index250coadd, ...
 coefs250coadd, Nalphacoadd, Nbetacoadd] = computeIndex(alpha, beta, ...
                                                        pointing250, ...
                                                        6, 6, ...
                                                        Nbolo250, ...
                                                        N_scan_total);

[alpha, beta, pointing250bis, index250, coefs250, Nalpha, Nbeta] = ...
    computeIndex(alpha, beta, pointing250, alpha_step, beta_step, ...
                 Nbolo250, N_scan_total);

%% The coaddition without offsets and estimation of offsets
offsets = zeros(1,Nbolo250);

coadd = dirtymap(data{1}, index250coadd, coefs250coadd, offsets, ...
                 Nalphacoadd, Nbetacoadd, N_scan_total, ...
                 unique_speed, the_speeds, Nbolo250);

[coadd2 offsets] = dirtymapOffsets(data{1}, index250coadd, ...
                                   coefs250coadd, Nalphacoadd, ...
                                   Nbetacoadd, N_scan_total, ...
                                   unique_speed, the_speeds, ...
                                   Nbolo250, 20);
              
%% Inversion
init = zeros(Nalpha, Nbeta);
objectMean = zeros(size(init));

[Hrond250 Hdirect250] = computeRI(10, 1, band_250, ...
                                  central_wavelength_250, ...
                                  sigma_coef, sigma_alpha, ...
                                  sigma_beta, alpha_step, beta_step, ...
                                  Nalpha, Nbeta, Nspeed, Norder, ...
                                  unique_speed, opt_efficiency, ...
                                  time_constante, gain);

break


%% Gain
G = Hrond250(1);

[diffAlpha diffBeta circDalpha circDbeta] = computeReg(sigma_alpha, ...
                                                       sigma_beta, ...
                                                       sky_alpha_period, ...
                                                       sky_beta_period, ...
                                                       5, Nalpha, ...
                                                       Nbeta);

regOp = circDalpha + circDbeta;
                  
%% GPAC Options
OPS(1)	 = 5;		% Affichage toutes les OPS(1) itérations
OPS(2)	 = 1e-8;	% Seuil d'arret sur x : 1e-3
OPS(3)	 = 1e-8;	% Seuil d'arret sur f : 1e-3
OPS(4)	 = 15;		% Un coup de grad tous les OPS(4)
OPS(5)	 = 3;		% Strategie 3 == gradient conjugué (0,1,2,3)
OPS(6)	 = 1;		% Norme (0,1,2)
OPS(7)	 = 0;		% Minimisation en ligne (0,1)
OPS(14)	 = 40;		% Nbre maximum d'itérations
OPS(15)	 = 1e-10;	% Pas minimum (MuEps)
OPS(18)	 = 1e-6;		% Pas a l'origine

%% Conjugate gradient options
cgoptions.thresold = 1e-1000; cgoptions.maxIter = OPS(14); %cgoptions.numfig = 1000;

%% Hypers parameters value
hypers = zeros(2,1);
hypers(1,1) = 1e-5;
hypers(2,1) = 1e1;

%% Inversion with GPAC
tic
[mapGpac ooptim] = invGpac(init, OPS, data{1}, hypers, Hrond250, ...
                           index250, zeros(1,Nbolo250), regOp, ...
                           objectMean, Nalpha, Nbeta, Norder, ...
                           N_scan_total, Nbolo250, Nspeed, ...
                           unique_speed, the_speeds, 'fig',1000);
TGpac = toc;

%% Avec restimation des offsets venant de la coadd
tic
[mapGpacO offsetsGpac ooptimGpac] = invGpacOffsets(init, OPS, ...
                                                  data{1}, hypers, ...
                                                  Hrond250, ...
                                                  index250, offsets, ...
                                                  regOp, objectMean, ...
                                                  Nalpha, Nbeta, ...
                                                  Norder, ...
                                                  N_scan_total, ...
                                                  Nbolo250, Nspeed, ...
                                                  unique_speed, ...
                                                  the_speeds, 5);
TGpacO = toc;

%% Avec le gradient conjugué

tic
[mapCg ooptimCg] = inversion(init, cgoptions, data{1}, hypers, ...
                             Hrond250, index250, coefs250, offsets, ...
                             regOp, objectMean, Nalpha, Nbeta, ...
                             Norder, N_scan_total, Nbolo250, Nspeed, ...
                             unique_speed, the_speeds);
TCg = toc;

%% Avec le gradient conjugué et estimation des offsets

tic
[mapCgO offsetsCg ooptimCg] = inversionOffsets(init, cgoptions, ...
                                              data{1}, hypers, ...
                                              Hrond250, index250, ...
                                              coefs250, offsetsCg, ...
                                              regOp, objectMean, ...
                                              Nalpha, Nbeta, Norder, ...
                                              N_scan_total, ...
                                              Nbolo250, Nspeed, ...
                                              unique_speed, the_speeds, 1);
TCgO = toc;

critGpac = calcQuadCrit(mapGpacO, data{1}, hypers, Hrond250, ...
                        index250, regOp, objectMean, Nalpha, Nbeta, ...
                        Norder, Nbolo250, Nspeed, N_scan_total, ...
                        unique_speed, the_speeds)
                           
critCg = calcQuadCrit(mapCgO, data{1}, hypers, Hrond250, index250, ...
                      regOp, objectMean, Nalpha, Nbeta, Norder, ...
                      Nbolo250, Nspeed, N_scan_total, unique_speed, ...
                      the_speeds)
                           
disp(['Speed gain : ',num2str(TGpac/TCg)])
disp(['Speed gain w o : ',num2str(TGpacO/TCgO)])

figure(1)
clf

subplot(121)
imagesc(alphaCoadd, betaCoadd, coadd2/G); axis image; colormap(hot); colorbar

subplot(122)
imagesc(alpha, beta, mapCgO); axis image; colormap(hot); colorbar

figure(2)
clf

subplot(121)
imagesc(alphaCoadd, betaCoadd, coadd2/G); axis image; colormap(cm); colorbar

subplot(122)
imagesc(alpha, beta, mapCgO); axis image; colormap(cm); colorbar

figure(3)
clf

subplot(121)
imagesc(alphaCoadd, betaCoadd, log(coadd2/G + abs(min(coadd2(:)/G)))); axis image; colormap(hot); colorbar

subplot(122)
imagesc(alpha, beta, log(mapCgO + abs(min(mapCgO(:))))); axis image; colormap(hot); colorbar

figure(4)
theAlpha = -100;
ligne = find(alphaCoadd <= theAlpha, 1, 'last' );
ligne2 = find(alpha <= theAlpha, 1, 'last' );
clf
plot(beta, mapCgO(ligne2,:),'r')
hold on
plot(betaCoadd, coadd2(ligne,:)/G)

break

%% write the fits file
path(path,'../mfitsio-1.2.4-src/mfitsio/')
map=fliplr(init(:,:,1));
%map=init(:,:,1);
path(path,'../')
make_header

fits_write(fits_name,S,C,map');
