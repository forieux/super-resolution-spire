set(0,'defaultaxesfontsize',25);
set(0,'defaulttextfontsize',25);
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);

SizeMarker = 2;
SizeLabel = 30;

figplace = './figs2/';

%%
load /mnt/espace/espace/results/aa2/doux_mean_above_true_s1e4.mat
sc = [min(skyEap(:)) max(skyEap(:))];

figure(1)
clf()
imagesc(skyEap, sc)
colormap(gray)
colorbar()
axis xy; axis square

print('-depsc',[figplace,'/doux_mean_above_true_s1e4']);

%%
figure(1)
clf()
plot(instChain,'.')
hold on
plot(ones(length(instChain),1)*true_sigma_coef, 'r--')
plot(ones(length(instChain),1)*meanInst, 'k--')
plot(ones(length(instChain),1)*(meanInst + sigmaInst), 'k.')
plot(ones(length(instChain),1)*(meanInst - sigmaInst), 'k.')
grid on
ylim([1.5e4, 6e4])
legend('Chain', 'True', 'Mean', 'Mean+-sigma')

print('-depsc',[figplace,'/doux_mean_above_true_s1e4_chain']);

%%
disp('doux_mean_above_true_s1e4')
disp(['Vrai :', num2str(true_sigma_coef)])
disp(['Mean (sigma) :', num2str(meanInst), ' (', num2str(sigmaInst),')'])
disp(['est - true :', num2str(mean(instChain(500:end)) - true_sigma_coef)])
clear std
disp(['Var :', num2str(std(instChain(500:end)))])

%%
load /mnt/espace/espace/results/aa2/doux_true_above_mean_s1e4.mat

figure(1)
clf()
imagesc(skyEap, sc)
colormap(gray)
colorbar()
axis xy; axis square

print('-depsc',[figplace,'/doux_true_above_mean_s1e4']);

%%
figure(1)
clf()
plot(instChain,'.')
hold on
plot(ones(length(instChain),1)*true_sigma_coef, 'r--')
plot(ones(length(instChain),1)*meanInst, 'k--')
plot(ones(length(instChain),1)*(meanInst + sigmaInst), 'k.')
plot(ones(length(instChain),1)*(meanInst - sigmaInst), 'k.')
grid on
ylim([1.5e4, 6e4])
legend('Chain', 'True', 'Mean', 'Mean+-sigma')

print('-depsc',[figplace,'/doux_true_above_mean_s1e4_chain']);

%%
disp('doux_true_above_mean_s1e4')
disp(['Vrai :', num2str(true_sigma_coef)])
disp(['Mean (sigma) :', num2str(meanInst), ' (', num2str(sigmaInst),')'])
disp(['est - true :', num2str(mean(instChain(500:end)) - true_sigma_coef)])
clear std
disp(['Var :', num2str(std(instChain(500:end)))])


%%
load /mnt/espace/espace/results/aa2/pics_mean_above_true_s1e4.mat

figure(1)
clf()
imagesc(skyEap, sc)
colormap(gray)
colorbar()
axis xy; axis square

print('-depsc',[figplace,'/pics_mean_above_true_s1e4']);

%%
figure(1)
clf()
plot(instChain,'.')
hold on
plot(ones(length(instChain),1)*true_sigma_coef, 'r--')
plot(ones(length(instChain),1)*meanInst, 'k--')
plot(ones(length(instChain),1)*(meanInst + sigmaInst), 'k.')
plot(ones(length(instChain),1)*(meanInst - sigmaInst), 'k.')
grid on
ylim([1.5e4, 6e4])
legend('Chain', 'True', 'Mean', 'Mean+-sigma')

print('-depsc',[figplace,'/pics_mean_above_true_s1e4_chain']);


%%
disp('pics_mean_above_true_s1e4')
disp(['Vrai :', num2str(true_sigma_coef)])
disp(['Mean (sigma) :', num2str(meanInst), ' (', num2str(sigmaInst),')'])
disp(['est - true :', num2str(mean(instChain(500:end)) - true_sigma_coef)])
disp(['Est :', num2str(mean(instChain(500:end)))])
clear std
disp(['Var :', num2str(std(instChain(500:end)))])


%%
load /mnt/espace/espace/results/aa2/pics_true_above_mean_s1e4.mat

figure(1)
clf()
imagesc(skyEap, sc)
colormap(gray)
colorbar()
axis xy; axis square

print('-depsc',[figplace,'/pics_true_above_mean_s1e4']);

%%
figure(1)
clf()
plot(instChain,'.')
hold on
plot(ones(length(instChain),1)*true_sigma_coef, 'r--')
plot(ones(length(instChain),1)*meanInst, 'k--')
plot(ones(length(instChain),1)*(meanInst + sigmaInst), 'k.')
plot(ones(length(instChain),1)*(meanInst - sigmaInst), 'k.')
grid on
ylim([1.5e4, 6e4])
legend('Chain', 'True', 'Mean', 'Mean+-sigma')

print('-depsc',[figplace,'/pics_true_above_mean_s1e4_chain']);


%%
disp('pics_true_above_mean_s1e4')
disp(['Vrai :', num2str(true_sigma_coef)])
disp(['Mean (sigma) :', num2str(meanInst), ' (', num2str(sigmaInst),')'])
disp(['est - true :', num2str(mean(instChain(500:end)) - true_sigma_coef)])
clear std
disp(['Var :', num2str(std(instChain(500:end)))])
