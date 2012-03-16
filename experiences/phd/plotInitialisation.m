figplace = '/home/orieux/tex/manuscrit/imagesRes';

set(0,'defaulttextinterpreter','none')

set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxeslinewidth',1);
set(0,'defaultlinelinewidth',1);

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
imagesc(SupAlpha, SupBeta, Hdirect360(:,:,1,1)); setim
ylabel('$\alpha$')
xlabel('$\beta$')

name = 'RI360ordre0';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
ligne = floor(Nsalpha/2);
plot(SupAlpha, Hdirect360(:,ligne,1,1))
xlabel('$\alpha$ (avec $\beta = 0$)')
grid on

name = 'sliceRI360ordre0';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
H = circshift(Hrond360(:,:,1,1), [floor(Nsalpha/2) floor(Nsbeta/2)]);
imagesc(FALPHAS(1:Nsalpha,1), FBETAS(1,1:Nsbeta), abs(H(1:Nsalpha, 1:Nsbeta)))
setim;
xlabel('$\fb$')
ylabel('$\fa$')

name = 'TF360ordre0';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
imagesc(SupAlpha, SupBeta, Hdirect360(:,:,2,1)); setim
ylabel('$\alpha$')
xlabel('$\beta$')

name = 'RI360ordre1';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
ligne = floor(Nsalpha/2);
plot(SupAlpha, Hdirect360(:,ligne,2,1))
xlabel('$\alpha$ (avec $\beta = 0$)')
grid on

name = 'sliceRI360ordre1';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
H = circshift(Hrond360(:,:,2,1), [floor(Nsalpha/2) floor(Nsbeta/2)]);
imagesc(FALPHAS(1:Nsalpha,1), FBETAS(1,1:Nsbeta), abs(H(1:Nsalpha, 1: Nsbeta)))
setim; 
xlabel('$\fb$')
ylabel('$\fa$')

name = 'TFRI360ordre1';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%% Pointing

color = ['b' 'g' 'r' 'c' 'm' 'y' 'k'];

clf
for iscan = 1:6
    p = oldPointing360{iscan};
    plot(p(1,1:10:end),p(2,1:10:end),['.',color(iscan)],'MarkerSize',SizeMarker);
    hold on
end
axis image

name = 'pointage';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
p = oldPointing360{1};
p2 = pointing360{1};

plot(p(1,:),p(2,:),'.','MarkerSize',SizeMarker);
hold on
plot(p2(1,:),p2(2,:),'.k','MarkerSize',SizeMarker);
axis image

xlim([580 660])
ylim([330 410])

xlabel('$\beta$')
ylabel('$\alpha$')

name = 'approxPointage';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
coefs = zeros(Nalpha, Nbeta);
for iscan = 1:6
    coefs = coefs + coefs250{iscan};
    imagesc(alpha, beta, coefs);
    setim
end
minCoef = min(coefs(:));
maxCoef = max(coefs(:));

xlabel('$\beta$')
ylabel('$\alpha$')

name = 'redondance250';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
coefs = zeros(Nalpha, Nbeta);
for iscan = 1:6
    coefs = coefs + coefs360{iscan};
    imagesc(alpha, beta, coefs, [minCoef maxCoef]);
    setim
end

xlabel('$\beta$')
ylabel('$\alpha$')

name = 'redondance360';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

clf
coefs = zeros(Nalpha, Nbeta);
for iscan = 1:6
    coefs = coefs + coefs520{iscan};
    imagesc(alpha, beta, coefs, [minCoef maxCoef]);
    setim
end

xlabel('$\beta$')
ylabel('$\alpha$')

name = 'redondance520';
laprint(1,[figplace,'/',name], 'width',8)
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
xlabel('$t$ [s]')
ylabel('$\yb$ [V]')

name = 'donnees';
laprint(1,[figplace,'/',name], 'width',10)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%

load ../data&true/simul_ism_SPS_comp
sky2 = zeros(Nalpha, Nbeta);
ma = min([Nalpha, size(sky,1)]);
mb = min([Nbeta,  size(sky,2)]);
sky2(1:ma,1:mb) = sky(1:ma,1:mb);
sky2 = fftshift(circshift(sky2,[-floor(ma/2) -floor(mb/2)]));
sky = sky2; clear sky2;
sky = sky*1e-4;

clf
imagesc(alpha, beta, sky); setim

name = 'cirrus';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%

load ../data&true/simul_cond_ism_SPS_comp
sky2 = zeros(Nalpha, Nbeta);
ma = min([Nalpha, size(sky,1)]);
mb = min([Nbeta,  size(sky,2)]);
sky2(1:ma,1:mb) = sky(1:ma,1:mb);
sky2 = fftshift(circshift(sky2,[-floor(ma/2) -floor(mb/2)]));
sky = sky2; clear sky2;
sky = sky*1e-4;

clf
imagesc(alpha, beta, sky); setim

name = 'cirrusDot';
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

clf
imagesc(alpha, beta, sky); setim

name = 'galaxie';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

clf
imagesc(alpha, beta, sky, [0 8e-5]); setim

name = 'galaxieSat';
laprint(1,[figplace,'/',name], 'width',8)
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%


close all
