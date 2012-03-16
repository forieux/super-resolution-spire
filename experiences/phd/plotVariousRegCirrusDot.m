%% variousRegCirrusDCT - Compute several estimation with different reg

%% Base
% format long
path(path,'../libspire')
path(path,'../simulator')
path(path,'../utils')
path(path,'../')
path(path,'/home/seismic/matlab/optimgpi')

placemount = '/espace/orieux/results/'
expname = 'variousRegCirrusDotDCT';

Ntasks = 144;

variousReg = logspace(7,15,Ntasks);

name = [placemount,expname,'/initialization'];
load(name, 'sky')

% %% Film
% for itask = 1:length(variousReg)
    
%     name = [placemount,expname,'/gpac_',num2str(itask)];
%     load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')
    
%     imagesc(xchap(:,:,2)); setim
    
%     drawnow
    
%     itask
%     pause(0.3)
    
% end

%% L2
for itask = 1:length(variousReg)
    
    name = [placemount,expname,'/gpac_',num2str(itask)];
    load(name, 'xchap', 'hypers', 'SortieOPS', 'Histo')
    
    dist(itask) = sum2((sky(:,:,2) - xchap(:,:,2))); setim
    regul(itask) = hypers(2,2,1);
 
end

% For sorting
sorted= [regul',dist'];
sorted = sortrows(sorted,1);

% semilogx(sorted(:,1),sorted(:,2))

sorted(find(sorted(:,2) == min(sorted(:,2))),1)
