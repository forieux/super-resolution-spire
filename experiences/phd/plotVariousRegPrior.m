%% Base
% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

%%

placemount = '/espace/orieux/results/'
expname = 'variousRegPriorDCT'

Ntasks = 144;
variousReg = logspace(7,15,Ntasks);

name = [placemount,expname,'/initialization'];
load(name, 'sky')

% for itask = 1:length(variousReg)
    
%     name = [placemount,expname,'/gpac_',num2str(itask)];
%     load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')
    
%     imagesc(xchap(:,:,2)); setim
    
%     drawnow
    
%     itask
%     pause(0.3)
    
% end

%% Gaussian decomposition function
alphaS = [-5:5];
[ALPHAS BETAS] = ndgrid(alphaS,alphaS);

sigma_alpha = 0.6*2;
sigma_beta  = 0.6*2;

%% Fonction de décomposition 
fgaussian = 1/(2*pi*sigma_alpha*sigma_beta)* exp(-ALPHAS.^2/(2*sigma_alpha^2) ...
                                                 -  BETAS.^2/(2* ...
                                                  sigma_beta^2));

%% L2
for itask = 1:length(variousReg)
    
    name = [placemount,expname,'/gpac_',num2str(itask)];
    load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')
    
    dist2(itask) = sum2((conv2(sky(:,:,2) - xchap(:,:,2),fgaussian,'same')).^2)/sum2(conv2(sky(:,:,2),fgaussian,'same').^2);
    dist1(itask) = sum2(abs(conv2(sky(:,:,2) - xchap(:,:,2),fgaussian,'same')))/sum2(abs(conv2(sky(:,:,2),fgaussian,'same')));
    regul(itask) = hypers(2,2,1);
    
end
% For sorting
sorted1 = [regul',dist1'];
sorted1 = sortrows(sorted1,1);

sorted2 = [regul',dist2'];
sorted2 = sortrows(sorted2,1);

% semilogx(sorted(:,1),sorted(:,2))

theBest1 = find(sorted1(:,2) == min(sorted1(:,2)));
theBest2 = find(sorted2(:,2) == min(sorted2(:,2)));

semilogx(sorted1(:,1),sorted1(:,2))
plot(sorted1(theBest1,1),sorted1(theBest1,2),'.')
hold on
semilogx(sorted2(:,1),sorted2(:,2),'r')
plot(sorted2(theBest2,1),sorted2(theBest2,2),'.r')
