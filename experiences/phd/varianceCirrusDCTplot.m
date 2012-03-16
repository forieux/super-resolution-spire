numfig = 0;

set(0, 'defaultaxesfontsize', 25);
set(0, 'defaulttextfontsize', 20);
set(0, 'defaultaxeslinewidth', 2);
set(0, 'defaultlinelinewidth', 2);
set(0, 'DefaultTextInterpreter', 'tex')

SizeMarker = 15;
SizeLabel = 40;

%% Mean computation
skyMean = zeros(Nalpha,Nbeta,3);

for itask = 1:Ntasks
    
    for iworker = 1:length(alltasks)
       
        name = [placemount,expname,'/sample_',num2str(itask),'_',num2str(iworker)];
        load(name,'skySample')
        
        skyMean = skyMean + skySample;
        
    end
    
end

skyMean = skyMean / (Ntasks*8);

%% Var computation
skyVar = zeros(Nalpha,Nbeta,3);

for itask = 1:Ntasks
    
    for iworker = 1:length(alltasks)
       
        name = [placemount,expname,'/sample_',num2str(itask),'_',num2str(iworker)];
        load(name,'skySample')
        
        skyVar = skyVar + skySample.^2;%(skySample - skyMean).^2;
        
    end
    
end

skyVar = skyVar / (Ntasks*8);

skyVar = skyVar - skyMean.^2;

skyStd = sqrt(skyVar);

%% End
sol = conv2(skyMean(:,:,2),fgaussian,'same');
true = conv2(sky(:,:,2),fgaussian,'same');
solP3 = conv2(skyMean(:,:,2) + 3*skyStd(:,:,2), fgaussian, 'same');
solM3 = conv2(skyMean(:,:,2) - 3*skyStd(:,:,2), fgaussian, 'same');

solStd = conv2(skyStd(:,:,2),fgaussian,'same');

%%
plot(sol(200,:))
hold on
plot(sol(200,:) - solStd(200,:),'k')
plot(sol(200,:) + solStd(200,:),'k')

