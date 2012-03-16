clear all

format('long');

load ../data&true/data_fits_image
load ../data&true/pos_fits_image




 %% for polaris
% d=data{1};
% data{1}=d(1:64);
% 
% p=pointing{1};
% pointing{1}=p(1:64);

%% set path for output

load mycm

tic
racine = '../spire/';
expname = 'ngcPMW/';
system(['mkdir -p ',racine,expname]);
name = '/spire';



%% set path, specially gpac
path(path,'../')
path(path,'../libspire')
path(path,'../utils')

scan_tot = size(data{1},2);%10;
%pointing250 = pointing{1};
pointing_band = pointing{2};
% pointing520 = pointing{3};

temp_ra=[];
temp_dec=[];

%% determining ra_mean and dec_mean
for iscan = 1:scan_tot
    p = pointing_band{iscan};
    temp_ra=[temp_ra deg2rad(p(1,:))];
    temp_dec=[temp_dec deg2rad(p(2,:))];
end

% break;

raz=mean(temp_ra);
decz=mean(temp_dec);

clear temp_ra temp_dec

% tan_ra=deg2rad(315.47608333333335);
% tan_dec=deg2rad(68.15972222222223);

%% M�me pixel avant et apr�s projection pour la coordonn�e de r�f�rence
%% (quelque soit le syt�me de coordonn�e)
% ptot1=[];
% ptot2=[];
% for iscan = 1:scan_tot
%   p = pointing250{iscan};
%   ptot1 = union(p(1,:),ptot1);
%   ptot2 = union(p(2,:),ptot2);
% end
% ptot = [ptot1;ptot2];
% ptot = deg2rad(ptot);
% [xi,eta] = convert_ra_dec_bis(ptot, tan_ra, tan_dec );

%% Projection
for iscan = 1:scan_tot
    p = pointing_band{iscan};
    p(2,:) = deg2rad(p(2,:));
    p(1,:) = deg2rad(p(1,:));
    [xi,eta] = convert_ra_dec_bis( p, raz, decz);
    z(1,:) = rad2arcsec(xi);
    z(2,:) = rad2arcsec(eta);
    pointing_band(iscan) = {z};
    clear z
end

clear xi eta

paramsInstrumentNGC7023
paramsObsNGC7023
paramsSkyNGC7023

bound = 4*sigma_coef*max(band_520) + 4*60; % arcsec
[alpha beta] = computeAxis(pointing_band, sky_alpha_period, sky_beta_period, ...
    bound, N_scan_total);

%% Precomputation of index and redondancy
[alphaCoadd, betaCoadd, pointing360coadd, index360coadd, ...
    coefs360coadd, Nalphacoadd, Nbetacoadd] = computeIndex(alpha, beta, ...
    pointing_band, ...
    12, 12, ... %% 6, 6
    Nbolo360, ...
    N_scan_total);

[alpha, beta, pointing360bis, index360, coefs360, Nalpha, Nbeta] = ...
    computeIndex(alpha, beta, pointing_band, alpha_step, beta_step, ...
    Nbolo360, N_scan_total);

%% The coaddition without offsets and estimation of offsets
offsets = zeros(1,Nbolo360);

% the_speeds=the_speeds(:,:,65:end);
% unique_speed=unique_speed(:,3:4);
% N_speed=2;

figure(100)
[coaddPMW offsetsPMW] = dirtymapOffsets(data{2}, index360coadd, ...
    coefs360coadd, Nalphacoadd, ...
    Nbetacoadd, N_scan_total, ...
    unique_speed, the_speeds, ...
    Nbolo360, 20); % 20



% coaddtest = dirtymap(data{2}, index360coadd, coefs360coadd, offsets, Nalphacoadd, Nbetacoadd, ...
%     N_scan_total, unique_speed, the_speeds, Nbolo360);
% 
% figure
% imagesc(coaddtest);
figure
imagesc(coaddPMW);
break;

%% Inversion
init = zeros(Nalpha, Nbeta);
objectMean = zeros(size(init));

[Hrond360 Hdirect360] = computeRI(10, 1, band_360, ...
    central_wavelength_360, ...
    sigma_coef, sigma_alpha, ...
    sigma_beta, alpha_step, beta_step, ...
    Nalpha, Nbeta, Nspeed, Norder, ...
    unique_speed, opt_efficiency, ...
    time_constante, gain);
%% Gain
G = Hrond360(1);

[diffAlpha diffBeta circDalpha circDbeta] = computeReg(sigma_alpha, ...
    sigma_beta, ...
    sky_alpha_period, ...
    sky_beta_period, ...
    5, Nalpha, ...
    Nbeta);

regOp = circDalpha + circDbeta;

%% Conjugate gradient options
cgoptions.thresold = 1e-10; cgoptions.maxIter = 40; %cgoptions.numfig = 1000;

% hyp_tab=[2e-5 3e-5 5e-5 6e-5 5e-5 7e-5 9e-5];
hyp_tab=[1e-6 2e-6 3e-6 5e-6];

for hyp_ind=1:size(hyp_tab,2)
    
    hyp_ind
%% Hypers parameters value
hypers = zeros(2,1);
hypers(1,1) = hyp_tab(hyp_ind); % 1e-5
hypers(2,1) = 1e1;

fits_name = strcat('result_PMW_6arc_19mai',num2str(hypers(1,1)),'.fits');

[mapCgO offsetsCg ooptimCg] = inversionOffsets(init, cgoptions, ...
    data{2}, hypers, ...
    Hrond360, index360, ...
    coefs360, offsetsPMW, ... 
    regOp, objectMean, ...
    Nalpha, Nbeta, Norder, ...
    N_scan_total, ...
    Nbolo360, Nspeed, ...
    unique_speed, the_speeds, 1);

tt=toc;
disp(['Speed : ',num2str(tt)])

% figure(1)
% clf
% 
% subplot(121)
% imagesc(alphaCoadd, betaCoadd, coaddPMW/G); axis image; colormap(hot); colorbar
% 
% subplot(122)
% imagesc(alpha, beta, mapCgO); axis image; colormap(hot); colorbar
% 
% figure(2)
% clf
% 
% subplot(121)
% imagesc(alphaCoadd, betaCoadd, coaddPMW/G); axis image; colormap(cm); colorbar
% 
% subplot(122)
% imagesc(alpha, beta, mapCgO); axis image; colormap(cm); colorbar
% 
% figure(3)
% clf
% 
% subplot(121)
% imagesc(alphaCoadd, betaCoadd, log(coaddPMW/G + abs(min(coaddPMW(:)/G)))); axis image; colormap(hot); colorbar
% 
% subplot(122)
% imagesc(alpha, beta, log(mapCgO + abs(min(mapCgO(:))))); axis image; colormap(hot); colorbar
% 
% figure(4)
% theAlpha = -100;
% ligne = find(alphaCoadd <= theAlpha, 1, 'last' );
% ligne2 = find(alpha <= theAlpha, 1, 'last' );
% clf
% plot(beta, mapCgO(ligne2,:),'r')
% hold on
% plot(betaCoadd, coaddPMW(ligne,:)/G)




%% write the fits file
path(path,'/home/mhusson/Matfiles/mfitsio-1.2.4-src/mfitsio/')
map=flipud(mapCgO);
% map=mapCgO;
%map=init(:,:,1);
path(path,'/home/mhusson/Matfiles/')
make_header

fits_write(fits_name,S,C,map);

end
