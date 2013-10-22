clear all

format('long');

%load ../../../Data/data_polaris_legs_ordonne
%load ../../../Data/pointing_polaris_legs_ordonne

% load ../data&true/data_fits_image
% load ../data&true/pos_fits_image

load ../data&true/data_ngc_notime
load ../data&true/pos_ngc_notime

%load ../../../Data/data_pola_final
%load ../../../Data/pointing_pola_final
% load ../../../Data/data_fits_pola_legs_PSW
% load ../../../Data/pointing_fits_pola_legs_PSW
% load data_part1_polaris
% load pointing_part1_polaris


%% for polaris
% d=data{1};
% data{1}=d(1:64);
%
% p=pointing{1};
% pointing{1}=p(1:64);

load mycm

tic
racine = '../spire/';
expname = 'ngc/';
system(['mkdir -p ',racine,expname]);
name = '/spire';

%% set path, specially gpac
path(path,'../')
path(path,'../libspire')
path(path,'../utils')
path(path,'../')

scan_tot = size(data{1},2);%10;
pointing_band = pointing{1};
% pointing360 = pointing{2};
% pointing520 = pointing{3};

temporal_sampling_periode = 1/18.666;

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

paramsInstrumentNGC7023 % to get Nbolo before the projection loop

vit_alpha_mean=zeros(1,iscan);
vit_beta_mean=zeros(1,iscan);
dir_alpha=zeros(1,iscan);
dir_beta=zeros(1,iscan);

%% Projection
for iscan = 1:scan_tot
    p = pointing_band{iscan};
    p(2,:) = deg2rad(p(2,:));
    p(1,:) = deg2rad(p(1,:));
    [xi,eta] = convert_ra_dec_bis( p, raz, decz );
    z(1,:) = rad2arcsec(xi);
    z(2,:) = rad2arcsec(eta);
    pointing_band(iscan) = {z};
    z1_reshape = reshape( z(1,:) ,size(z(1,:),2)/Nbolo250,Nbolo250);
    z2_reshape = reshape( z(2,:) ,size(z(2,:),2)/Nbolo250,Nbolo250);

    zdiff1= z1_reshape(2:end,:) - z1_reshape(1:end-1,:);%diff(z1_reshape,1,1);
    zdiff2= z2_reshape(2:end,:) - z2_reshape(1:end-1,:);%diff(z2_reshape,1,1);
    vit_alpha=zdiff1./temporal_sampling_periode;
    vit_beta=zdiff2./temporal_sampling_periode;
    vit_alpha_mean(iscan)=mean(mean(vit_alpha));
    vit_beta_mean(iscan)=mean(mean(vit_beta));

    clear z zdiff1 zdiff2 vit_alpha vit_beta
end

clear xi eta

dir_alpha=sign(vit_alpha_mean);
dir_beta=sign(vit_beta_mean);


paramsObsNGC7023_real_speeds
paramsSkyNGC7023

bound = 4*sigma_coef*max(band_520) + 4*60; % arcsec
[alpha beta] = computeAxis(pointing_band, sky_alpha_period, sky_beta_period, ...
    bound, N_scan_total);

%% Precomputation of index and redondancy
[alphaCoadd, betaCoadd, pointing250coadd, index250coadd, ...
    coefs250coadd, Nalphacoadd, Nbetacoadd] = computeIndex(alpha, beta, ...
    pointing_band, ...
    12, 12, ... %% 6, 6
    Nbolo250, ...
    N_scan_total);

[alpha, beta, pointing250bis, index250, coefs250, Nalpha, Nbeta] = ...
    computeIndex(alpha, beta, pointing_band, alpha_step, beta_step, ...
    Nbolo250, N_scan_total);

%% The coaddition without offsets and estimation of offsets
offsets = zeros(1,Nbolo250);



[coaddPSW offsetsPSW] = dirtymapOffsets(data{1}, zeros(Nbolo250, 1), index250coadd, ...
    Nalphacoadd, ...
    Nbetacoadd, N_scan_total, ...
    Nbolo250, 20); % 20


% coaddtest = dirtymap(data{1}, index250coadd, offsetscoefs250coadd, offsets, Nalphacoadd, Nbetacoadd, ...
%     N_scan_total, unique_speed, the_speeds, Nbolo250);

figure(1)
imagesc(coaddPSW); axis image

break


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
%% Gain
G = Hrond250(1);

[diffAlpha diffBeta circDalpha circDbeta] = computeReg(sigma_alpha, ...
    sigma_beta, ...
    sky_alpha_period, ...
    sky_beta_period, ...
    5, Nalpha, ...
    Nbeta);

regOp = circDalpha + circDbeta;

%% Conjugate gradient options
cgoptions.thresold = 1e-10; cgoptions.maxIter = 40; %cgoptions.numfig = 1000;

%hyp_tab=[1e-7 2e-7 3e-7 4e-7 5e-7 6e-7 7e-7 8e-7 9e-7];
hyp_tab=[5e-7];

figure(2)
subplot(2,2,1)
imagesc(coefs250{1}); axis image

subplot(2,2,2)
imagesc(Hdirect250(:,:,1,3)); axis image

subplot(2,2,3)
imagesc(coefs250{6}); axis image

subplot(2,2,4)
imagesc(Hdirect250(:,:,1,4)); axis image

figure(3)
subplot(2,2,1)
imagesc(coefs250{16}); axis image

subplot(2,2,2)
imagesc(Hdirect250(:,:,1,1)); axis image

subplot(2,2,3)
imagesc(coefs250{15}); axis image

subplot(2,2,4)
imagesc(Hdirect250(:,:,1,2)); axis image

%break

data250=data{1};

for hyp_ind=1:size(hyp_tab,2)

    disp(['Hyper ',num2str(hyp_ind),'/',num2str(size(hyp_tab,2))])
    %% Hypers parameters value
    hypers = zeros(2,1);
    hypers(1,1) = hyp_tab(hyp_ind); % 1e-5
    hypers(2,1) = 1e1;

    fits_name = strcat('result_PSW_6arc_',num2str(10/hypers(1,1)),'_0_filter.fits');

    for iter=1:3

        iter
        %%% Estimation de x sachant Offset et y
        %     NewData=data;
        %
        %     for iscan = 1:10
        %         Un = ones(size(data250{iscan},1),1);
        %         a = data250{iscan} - Un*Offset250.';
        %         NewData250(iscan) = {a};
        %     end
        %     NewData(1)={NewData250};


        [mapCgO offsetsCg ooptimCg] = inversionOffsets(init, cgoptions, ...
            data{1}, hypers, ...
            Hrond250, index250, ...
            coefs250, offsets, ... %offsetsPSW
            regOp, objectMean, ...
            Nalpha, Nbeta, Norder, ...
            N_scan_total, ...
            Nbolo250, Nspeed, ...
            unique_speed, the_speeds, 1);


        SortieModel = directInvariant(mapCgO, Hrond250, index250, ...
            Nalpha, Nbeta, Norder, Nbolo250, Nspeed, scan_tot, unique_speed, the_speeds);



        Sortie250=  SortieModel;
        a = data250{1}-Sortie250{1};
        for iscan = 2:10
            b = data250{iscan}-Sortie250{iscan};
            a = cat(1,a,b);
        end
        offsets = mean(a,1);

    end
    tt=toc;
    disp(['Speed : ',num2str(tt)])

    % figure(1)
    % clf
    %
    % subplot(121)
    % imagesc(alphaCoadd, betaCoadd, coaddPSW/G); axis image; colormap(hot); colorbar
    %
    % subplot(122)
    % imagesc(alpha, beta, mapCgO); axis image; colormap(hot); colorbar
    %
    % figure(2)
    % clf
    %
    % subplot(121)
    % imagesc(alphaCoadd, betaCoadd, coaddPSW/G); axis image; colormap(cm); colorbar
    %
    % subplot(122)
    % imagesc(alpha, beta, mapCgO); axis image; colormap(cm); colorbar
    %
    % figure(3)
    % clf
    %
    % subplot(121)
    % imagesc(alphaCoadd, betaCoadd, log(coaddPSW/G + abs(min(coaddPSW(:)/G)))); axis image; colormap(hot); colorbar
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
    % plot(betaCoadd, coaddPSW(ligne,:)/G)

end
    %% write the fits file
        path(path,'/home/mhusson/Matfiles/mfitsio-1.2.4-src/mfitsio/')
        map=flipud(mapCgO);
        % map=mapCgO;
        %map=init(:,:,1);
        path(path,'/home/mhusson/Matfiles/')
        make_header

        fits_write(fits_name,S,C,map);


figure(5)
imagesc(mapCgO); axis image
