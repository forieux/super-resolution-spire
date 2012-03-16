%% Base
% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

load /espace/orieux/results/init.mat

figplace = '/home/orieux/tex/manuscrit/imagesRes';

set(0,'defaulttextinterpreter','none')

% set(0,'defaultaxesfontsize',18);
% set(0,'defaulttextfontsize',15);
% set(0,'defaultaxeslinewidth',1);
% set(0,'defaultlinelinewidth',1);

SizeMarker = 2;
SizeLabel = 30;

numfig = 1;
figure(numfig)

falpha = linspace(-0.5,0.5,Nalpha);
fbeta = linspace(-0.5,0.5,Nbeta);

[FALPHAS FBETAS] = ndgrid(falpha,fbeta);
FALPHAS = circshift(ifftshift(FALPHAS), [floor(Nsalpha/2) 1]);
FBETAS = circshift(ifftshift(FBETAS), [1 floor(Nsbeta/2)]);

%%
%%
%% RI

clf
subplot(2,1,1)
imagesc(SupAlpha, SupBeta, Hdirect360(:,:,1,3)); 
axis image
axis off
%axis xy
colormap(gray)

subplot(2,1,2)
ligne = floor(Nsalpha/2);
plot(SupBeta, Hdirect360(ligne,:,1,3))
axis square
grid on

name = 'RI360ordre0';
laprint(1,[figplace,'/',name], 'options','factory','width',4,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
H = circshift(Hrond360(:,:,1,3), [floor(Nsalpha/2) floor(Nsbeta/2)]);
subplot(2,1,1)
imagesc(FALPHAS(1:Nsalpha,1)/alpha_step, FBETAS(1,1:Nsbeta)/beta_step, abs(H(1:Nsalpha, 1:Nsbeta)))
axis image
axis off
axis xy
colormap(gray)

subplot(2,1,2)
plot(FALPHAS(1:Nsalpha,1)/alpha_step, abs(H(1:Nsalpha, floor(Nsbeta/2))))
axis square
grid on

name = 'TF360ordre0';
laprint(1,[figplace,'/',name], 'options','factory','width',4,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
subplot(2,1,1)
imagesc(SupAlpha, SupBeta, Hdirect360(:,:,2,3))
axis image
axis off
axis xy
colormap(gray)

ligne = floor(Nsalpha/2);
subplot(2,1,2)
plot(SupBeta, Hdirect360(ligne,:,2,3))
axis square
grid on

name = 'RI360ordre1';
laprint(1,[figplace,'/',name], 'options','factory','width',4,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
H = circshift(Hrond360(:,:,2,3), [floor(Nsalpha/2) floor(Nsbeta/2)]);
subplot(2,1,1)
imagesc(FALPHAS(1:Nsalpha,1)/alpha_step, FBETAS(1,1:Nsbeta)/beta_step, abs(H(1:Nsalpha, 1: Nsbeta)))
axis image
axis off
axis xy
colormap(gray)
 
subplot(2,1,2)
plot(FALPHAS(1:Nsalpha,1)/alpha_step, abs(H(1:Nsalpha, floor(Nsbeta/2))))
axis square
grid on

name = 'TFRI360ordre1';
laprint(1,[figplace,'/',name], 'options','factory','width',4,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%% Pointing

color = ['b' 'g' 'r' 'c' 'm' 'y' 'k'];

clf
subplot(1,2,1)
for iscan = 1:6
    p = oldPointing360{iscan};
    plot(p(1,1:10:end),p(2,1:10:end),['.',color(iscan)],'MarkerSize',SizeMarker);
    hold on
end
for iscan = 1:6
    h = arrow([start_end_position(1,1,iscan) start_end_position(2,1,iscan)], ...
              [start_end_position(1,2,iscan) start_end_position(2,2, iscan)]);
    arrow(h,'Length',20)
    hold on
end
axis image

ylabel('$\beta$')
xlabel('$\alpha$')

% name = 'pointage';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

p = oldPointing360{1};
p2 = pointing360{1};

subplot(1,2,2)
plot(p(1,:),p(2,:),'o','MarkerSize',4);
hold on
plot(p2(1,:),p2(2,:),'xk','MarkerSize',4);
grid on
axis image

xlim([635 645])
ylim([375 385])

ylabel('$\beta$')
xlabel('$\alpha$')

name = 'approxPointage';
laprint(1,[figplace,'/',name], 'options','factory','width',12,'factor',0.8,'scalefonts','on')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
coefs = zeros(Nalpha, Nbeta);
for iscan = 1:6
    coefs = coefs + coefs250{iscan};
end
coefs(1,:) = 6;
coefs(end,:) = 6;
coefs(:,1) = 6;
coefs(:,end) = 6;

imagesc(alpha, beta, coefs);
axis image
axis xy; axis off
colorbar
colormap(flipud(gray))

title('$250~\mu\textrm{m}$')

minCoef = min(coefs(:));
maxCoef = max(coefs(:));

ylabel('$\beta$')
xlabel('$\alpha$')

name = 'redondance250';
laprint(1,[figplace,'/',name], 'width',5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

coefs = zeros(Nalpha, Nbeta);
for iscan = 1:6
    coefs = coefs + coefs360{iscan};
end
coefs(1,:) = 6;
coefs(end,:) = 6;
coefs(:,1) = 6;
coefs(:,end) = 6;

imagesc(alpha, beta, coefs, [minCoef maxCoef]);
axis image
axis xy; axis off
colorbar
colormap(flipud(gray))

title('$360~\mu\textrm{m}$')

xlabel('$\alpha$')
ylabel('$\beta$')

name = 'redondance360';
laprint(1,[figplace,'/',name], 'width',5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

coefs = zeros(Nalpha, Nbeta);
for iscan = 1:6
    coefs = coefs + coefs520{iscan};
end

coefs(1,:) = 6;
coefs(end,:) = 6;
coefs(:,1) = 6;
coefs(:,end) = 6;

imagesc(alpha, beta, coefs, [minCoef maxCoef]);
axis image
axis xy; axis off
colorbar
colormap(flipud(gray))

title('$520~\mu\textrm{m}$')

xlabel('$\alpha$')
ylabel('$\beta$')

name = 'redondance520';
laprint(1,[figplace,'/',name], 'width',5)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

out360 = output{2};
out = out360{1};

noisy = data360{1};

time = [1:size(noisy,1)]*temporal_sampling_periode;

clf
plot(time,out(:,1:22:end))
hold on
plot(time,noisy(:,1:22:end))
grid on
xlim([min(time) max(time)])
ylim([0.03 0.065])
xlabel('$t$ [s]')
ylabel('$\yb$ [V]')

name = 'donnees';
laprint(1,[figplace,'/',name], 'options','factory','width',6,'factor',0.7,'scalefonts','on')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

name = [placemount,'variousRegCirrusDCT','/gpac_',num2str(1)];
load(name, 'xchap')
dynamique = [min2(xchap(:,:,2)) max2(xchap(:,:,2))];

load ../data&true/simul_ism_SPS_comp
sky2 = zeros(Nalpha, Nbeta);
ma = min([Nalpha, size(sky,1)]);
mb = min([Nbeta,  size(sky,2)]);
sky2(1:ma,1:mb) = sky(1:ma,1:mb);
sky2 = fftshift(circshift(sky2,[-floor(ma/2) -floor(mb/2)]));
sky = sky2; clear sky2;
sky = sky*1e-4;

clf
imagesc(alpha, beta, sky, dynamique);
axis image
axis xy
colorbar
colormap(gray)

name = 'cirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

minCirrus = min2(sky);
maxCirrus = max2(sky);

load ../data&true/simul_cond_ism_SPS_comp
sky2 = zeros(Nalpha, Nbeta);
ma = min([Nalpha, size(sky,1)]);
mb = min([Nbeta,  size(sky,2)]);
sky2(1:ma,1:mb) = sky(1:ma,1:mb);
sky2 = fftshift(circshift(sky2,[-floor(ma/2) -floor(mb/2)]));
sky = sky2; clear sky2;
sky = sky*1e-4;

% clf
% imagesc(alpha, beta, sky); 
% axis image
% axis xy
% colorbar
% colormap(gray)
 
% name = 'cirrusDot';
% laprint(1,[figplace,'/',name], 'width',8)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

clf
imagesc(alpha, beta, sky, dynamique); 
axis image
axis xy
colorbar
colormap(gray)
 
name = 'cirrusDotSat';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

load ../data&true/galaxie
sky2 = zeros(Nalpha, Nbeta);
ma = min([Nalpha, size(sky,1)]);
mb = min([Nbeta,  size(sky,2)]);
sky2(1:ma,1:mb) = sky(1:ma,1:mb);
sky2 = fftshift(circshift(sky2,[-floor(ma/2) -floor(mb/2)]));
sky = sky2; clear sky2;
sky = sky*1e-4;

% clf
% imagesc(alpha, beta, sky); 
% axis image
% axis xy
% colorbar
% colormap(gray)

% name = 'galaxie';
% laprint(1,[figplace,'/',name], 'width',6)
% movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

clf
imagesc(alpha, beta, sky, [0 8e-5]); 
axis image
axis xy
colorbar
colormap(gray)

name = 'galaxieSat';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%% Prior

gammaX = 4e11;
regOp = circDalpha + circDbeta;
%%% Same seed (see >> doc randn)
RandStream.setDefaultStream(RandStream('mt19937ar','seed',0));

sky = uifft2(ufft2(randn(Nalpha,Nbeta))./(sqrt(gammaX*regOp)));
%%%% positivity
sky = sky + abs(min2(sky));

sky = conv2(sky,fgaussian,'same');

clf
imagesc(alpha, beta, sky); 
axis image
axis xy
colorbar
colormap(gray)

name = 'imagePrior';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

close all
