set(0,'defaultaxesfontsize',25);
set(0,'defaulttextfontsize',25);
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);

SizeMarker = 2;
SizeLabel = 30;

figplace = '/home/orieux/tex/papers/aa/tex/figs2';

%%
load /mnt/space/results/aa2/doux_mean_above_true_s1e4.mat
sc = [min(skyEap(:)) max(skyEap(:))];

disp('doux_mean_above_true_s1e4')
disp(['Vrai :', num2str(true_sigma_coef)])
disp(['Mean (sigma) :', num2str(meanInst), ' (', num2str(sigmaInst),')'])
disp(['est - true :', num2str(mean(instChain(1000:end)) - true_sigma_coef)])
disp(['relative error :', num2str((mean(instChain(1000:end)) - true_sigma_coef) / true_sigma_coef)])
clear std
disp(['Var :', num2str(std(instChain(1000:end)))])
disp('------------------------------')

%%
load /mnt/space/results/aa2/doux_true_above_mean_s1e4.mat

disp('doux_true_above_mean_s1e4')
disp(['Vrai :', num2str(true_sigma_coef)])
disp(['Mean (sigma) :', num2str(meanInst), ' (', num2str(sigmaInst),')'])
disp(['est - true :', num2str(mean(instChain(1000:end)) - true_sigma_coef)])
disp(['relative error :', num2str((mean(instChain(1000:end)) - true_sigma_coef) / true_sigma_coef)])
clear std
disp(['Var :', num2str(std(instChain(1000:end)))])
disp('------------------------------')

%%
load /mnt/space/results/aa2/pics_mean_above_true_s1e4.mat

figure(1)
clf()
imagesc(skyEap, sc)
colormap(gray)
colorbar()
axis xy; axis square
print('-depsc',[figplace,'/pics_mean_above_true_s1e4']);

figure(1)
clf()
plot(instChain,'.')
hold on
plot(ones(length(instChain),1) * true_sigma_coef, 'r--')
plot(ones(length(instChain),1) * meanInst, 'k--')
plot(ones(length(instChain),1) * (meanInst + sigmaInst), 'k.')
plot(ones(length(instChain),1) * (meanInst - sigmaInst), 'k.')
grid on
ylim([1.5e4, 6e4])
legend('Chain', 'True', 'Mean', 'Mean+-sigma')
print('-depsc',[figplace,'/pics_mean_above_true_s1e4_chain']);

disp('pics_mean_above_true_s1e4')
disp(['Vrai :', num2str(true_sigma_coef)])
disp(['Mean (sigma) :', num2str(meanInst), ' (', num2str(sigmaInst),')'])
disp(['Est :', num2str(mean(instChain(1000:end)))])
disp(['est - true :', num2str(mean(instChain(1000:end)) - true_sigma_coef)])
disp(['relative error :', num2str((mean(instChain(1000:end)) - true_sigma_coef) / true_sigma_coef)])
clear std
disp(['Var :', num2str(std(instChain(1000:end)))])
disp('------------------------------')

%%
load /mnt/space/results/aa2/pics_true_above_mean_s1e4.mat

figure(1)
clf()
imagesc(skyEap, sc)
colormap(gray)
colorbar()
axis xy; axis square
print('-depsc',[figplace,'/pics_true_above_mean_s1e4']);

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

disp('pics_true_above_mean_s1e4')
disp(['Vrai :', num2str(true_sigma_coef)])
disp(['Mean (sigma) :', num2str(meanInst), ' (', num2str(sigmaInst),')'])
disp(['Est :', num2str(mean(instChain(1000:end)))])
disp(['est - true :', num2str(mean(instChain(1000:end)) - true_sigma_coef)])
disp(['relative error :', num2str((mean(instChain(1000:end)) - true_sigma_coef) / true_sigma_coef)])
clear std
disp(['Var :', num2str(std(instChain(1000:end)))])

%%
load ../ngc7023_unsupervised_with_offset_estimation.mat

figure(1)
clf
imagesc(skyEap)
axis image; colormap(hot); colorbar

clf

subplot(1,2,2, 'Position', [0.65, 0.1, 0.25, .85])
hist(offsetsChain(2,:), linspace(-20, 20, 40))
set(gca,'CameraUpVector', [1,0,0]);
xlim([-20 20])
set(gca, 'ytick', []);
grid on
subplot(1,2,1, 'Position', [0.1, 0.1, 0.6, .85])
plot(offsetsChain(2,:), '.')
hold on
plot(cumsum(offsetsChain(2,:)) ./ cumsum(ones(size(offsetsChain(2,:)))), 'r')
ylim([-20 20])
grid on
print('-depsc',[figplace,'/offset2']);

clf
subplot(1,2,2, 'Position', [0.65, 0.1, 0.25, .85])
hist(offsetsChain(54,:), linspace(-20, 20, 40))
set(gca,'CameraUpVector', [1,0,0]);
xlim([-20 20])
set(gca, 'ytick', []);
grid on
subplot(1,2,1, 'Position', [0.1, 0.1, 0.6, .85])
plot(offsetsChain(54,:), '.')
hold on
plot(cumsum(offsetsChain(54,:)) ./ cumsum(ones(size(offsetsChain(54,:)))), 'r')
ylim([-20 20])
grid on
print('-depsc',[figplace,'/offset54']);

%%

clf
semilogy(squeeze(gxChain(:)))
hold on
%ylim([min(gxChain(:))-1e11 max(gxChain(:))+1e11])
grid on
print('-depsc',[figplace,'/ngc-gxchain-jfg'])

clf
plot(squeeze(gnChain(:)))
hold on
ylim([0 0.2])
grid on
print('-depsc',[figplace,'/ngc-gnchain-jfg'])
