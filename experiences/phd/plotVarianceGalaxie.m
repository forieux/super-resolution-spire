%% variousRegCirrusDCT - Compute several estimation with different reg

%% Base
% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

%%

placemount = '/espace/orieux/results/'
expname = 'varianceGalaxieDCT';

Ntasks = 35;

load([placemount,expname,'/initialization']);

%% Mean

mean = zeros(Nalpha, Nbeta);

for itask = 1:Ntasks
    
    for iworker = 1:8
        name = [placemount,expname,'/sample_',num2str(itask),'_',num2str(iworker)];
        load(name, 'skySample')
        
        mean = mean + conv2(skySample(:,:,2),fgaussian,'same');
    
    end

end

mean = mean/(Ntasks*8);

%% Var

cumulant2 = zeros(Nalpha, Nbeta);

for itask = 1:Ntasks
    
    for iworker = 1:8
        name = [placemount,expname,'/sample_',num2str(itask),'_',num2str(iworker)];
        load(name, 'skySample')
        
        cumulant2 = cumulant2 + (conv2(skySample(:,:,2),fgaussian,'same')).^2;
    
    end

end
cumulant2 = cumulant2/(Ntasks*8);

variance = cumulant2 - mean.^2;
standartDev = sqrt(variance);

%% Plot

plot(mean(200,:))
plot(mean(200,280:300) + 3*variance(200,280:300),'k')
plot(mean(200,280:300) - 3*variance(200,280:300),'k')



