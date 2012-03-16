clear all

% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

figplace = '/home/orieux/tex/manuscrit/imagesRes/cirrus/';
system(['mkdir -p ',figplace,'/eps']);

set(0,'defaulttextinterpreter','none')

set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxeslinewidth',1);
set(0,'defaultlinelinewidth',0.5);

SizeMarker = 2;
SizeLabel = 30;

numfig = 1;
figure(numfig)

%%
%% Biais

boundAlpha = [160 420];
boundBeta = [160 420];

boundBetal = [160 420];

% disp('Biais')

Ntasks = 104;

cumul1 = zeros(Nalpha, Nbeta);
cumul2 = zeros(Nalpha, Nbeta);

for itask = 1:Ntasks
    name = [placemount,'biaisCirrusDCT','/gpac_',num2str(itask)];
    load(name, 'xchap')
    
    cumul1 = cumul1 + xchap(:,:,2);
    cumul2 = cumul2 + xchap(:,:,2).^2;
end
cumul1 = cumul1/Ntasks;
cumul2 = cumul2/Ntasks;

biais = cumul1 - sky(:,:,2);
variance = cumul2 - cumul1.^2;

clf

minBiais = -7e-6;
maxBiais = 7e-6;
scBiais = [minBiais maxBiais];

scVar = [0 4.2e-13];
scVar3 = [0 4.2e-11];
scVar2 = [0 3.2e-15];

subplot(211)
imagesc(alpha(boundAlpha(1):boundAlpha(2)),beta(boundBeta(1): boundBeta(2)), ...
        biais(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)), scBiais)
axis image
axis xy; axis off
colormap(gray)

subplot(212)
plot(beta(boundBetal(1):boundBetal(2)), biais(ligne, boundBetal(1): ...
                                              boundBetal(2)))
ylim(scBiais)
xlim([beta(boundBetal(1)) beta(boundBetal(2))])
axis square

grid on

name = 'biais360Cirrus';
laprint(1,[figplace,'/',name], 'options','factory','width',5,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% Variance

subplot(211)
imagesc(alpha(boundAlpha(1):boundAlpha(2)),beta(boundBeta(1): boundBeta(2)), ...
        variance(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)), scVar)
axis image
axis xy; axis off
colormap(gray)

subplot(212)
plot(beta(boundBetal(1): boundBetal(2)), variance(ligne, boundBetal(1): ...
                                                  boundBetal(2)))
ylim(scVar)
xlim([beta(boundBetal(1)) beta(boundBetal(2))])
axis square
grid on

name = 'varEstimateur360Cirrus';
laprint(1,[figplace,'/',name], 'options','factory','width',5,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

% 
% 
% Biais low reg

disp('Biais')

Ntasks = 104;

cumul1 = zeros(Nalpha, Nbeta);
cumul2 = zeros(Nalpha, Nbeta);

for itask = 1:Ntasks
    name = [placemount,'biaisCirrusDCTlowReg','/gpac_',num2str(itask)];
    load(name, 'xchap')
    
    cumul1 = cumul1 + xchap(:,:,2);
    cumul2 = cumul2 + xchap(:,:,2).^2;
end
cumul1 = cumul1/Ntasks;
cumul2 = cumul2/Ntasks;

biais = cumul1 - sky(:,:,2);
variance = cumul2 - cumul1.^2;

clf

subplot(211)
imagesc(alpha(boundAlpha(1):boundAlpha(2)),beta(boundBeta(1): boundBeta(2)), ...
        biais(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)), scBiais)
axis image
axis xy; axis off
colormap(gray)

subplot(212)
plot(beta(boundBetal(1): boundBetal(2)), biais(ligne, boundBetal(1):boundBetal(2)))
ylim(scBiais)
xlim([beta(boundBetal(1)) beta(boundBetal(2))])
axis square
grid on

name = 'lowRegBiais360Cirrus';
laprint(1,[figplace,'/',name], 'options','factory','width',5,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% Variance

subplot(211)
imagesc(alpha(boundAlpha(1):boundAlpha(2)),beta(boundBeta(1): boundBeta(2)), ...
        variance(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)), scVar3)
axis image
axis xy; axis off
colormap(gray)

subplot(212)
plot(beta(boundBetal(1): boundBetal(2)), variance(ligne, boundBetal(1):boundBetal(2)))
ylim(scVar3)
xlim([beta(boundBetal(1)) beta(boundBetal(2))])
axis square
grid on

name = 'lowRegvarEstimateur360Cirrus';
laprint(1,[figplace,'/',name], 'options','factory','width',5,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%%
%% Biais high reg 

disp('Biais')

Ntasks = 104;

cumul1 = zeros(Nalpha, Nbeta);
cumul2 = zeros(Nalpha, Nbeta);

for itask = 1:Ntasks
    name = [placemount,'biaisCirrusDCThighReg','/gpac_',num2str(itask)];
    load(name, 'xchap')
    
    cumul1 = cumul1 + xchap(:,:,2);
    cumul2 = cumul2 + xchap(:,:,2).^2;
end
cumul1 = cumul1/Ntasks;
cumul2 = cumul2/Ntasks;

biais = cumul1 - sky(:,:,2);
variance = cumul2 - cumul1.^2;

clf

subplot(211)
imagesc(alpha(boundAlpha(1):boundAlpha(2)),beta(boundBeta(1): boundBeta(2)), ...
        biais(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)), scBiais)
axis image
axis xy; axis off
colormap(gray)

subplot(212)
plot(beta(boundBetal(1): boundBetal(2)), biais(ligne, boundBetal(1): ...
                                               boundBetal(2)))
ylim(scBiais)
xlim([beta(boundBetal(1)) beta(boundBetal(2))])
axis square
grid on

name = 'highRegBiais360Cirrus';
laprint(1,[figplace,'/',name], 'options','factory','width',5,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])

%%
%% Variance

subplot(211)
imagesc(alpha(boundAlpha(1):boundAlpha(2)),beta(boundBeta(1): boundBeta(2)), ...
        variance(boundAlpha(1):boundAlpha(2),boundBeta(1):boundBeta(2)), scVar2)
axis image
axis xy; axis off
colormap(gray)

subplot(212)
plot(beta(boundBetal(1): boundBetal(2)), variance(ligne, boundBetal(1): ...
                                                  boundBetal(2)))
ylim(scVar2)
xlim([beta(boundBetal(1)) beta(boundBetal(2))])
axis square
grid on

name = 'highRegvarEstimateur360Cirrus';
laprint(1,[figplace,'/',name], 'options','factory','width',5,'factor',0.4,'scalefonts','off')
movefile([figplace,'/',name,'.eps'],[figplace,'/eps/',name,'.eps'])


close all

