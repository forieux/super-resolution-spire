path= '/Users/orieux/tmp/';

set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxeslinewidth',1);
set(0,'defaultlinelinewidth',1);

base = 'ngc7023';

sc = [min(mapCgO(:)) max(mapCgO(:))];

figure(1)
% clf
% imagesc(alphaCoadd, betaCoadd, coaddPSW/G)
% axis image; axis off; colormap(hot); 
% print('-depsc',[path,base,'-coadd']);
% 
% clf
% imagesc(alpha, beta, mapCgO)
% axis image; axis off; colormap(hot); 
% print('-depsc',[path,base,'-prop']);

% clf
% imagesc(alphaCoadd, betaCoadd, coaddPSW/G)
% axis image; axis off; colormap(cm); 
% print('-depsc',[path,base,'-coadd-sat']);
% 
% clf
% imagesc(alpha, beta, mapCgO)
% axis image; axis off; colormap(cm);
% print('-depsc',[path,base,'-prop-sat']);

clf
imagesc(alphaCoadd, betaCoadd, coaddPSW/G, sc)
axis image; axis off; colormap(cm); 
xlim([-400 400])
ylim([-400 400])
print('-depsc',[path,base,'-zoom-coadd-sat']);

clf
imagesc(alpha, beta, mapCgO, sc)
axis image; axis off; colormap(cm);
xlim([-400 400])
ylim([-400 400])
print('-depsc',[path,base,'-zoom-prop-sat']);

clf
% imagesc(alphaCoadd, betaCoadd, sqrt(coaddPSW/G + abs(min(mapCgO(:)))), sqrt(sc))
% axis image; axis off; colormap(hot); 
% xlim([-800 800])
% ylim([-800 800])

clf
imagesc(alphaCoadd, betaCoadd, coaddPSW/G, sc)
axis image; axis off; colormap(hsv); 
xlim([-800 800])
ylim([-800 800])
print('-depsc',[path,base,'-coadd-sq']);

clf
imagesc(alpha, beta, sqrt(mapCgO + abs(min(mapCgO(:)))), sqrt(sc))
axis image; axis off; colormap(hot);
xlim([-800 800])
ylim([-800 800])

clf
imagesc(alphaCoadd, betaCoadd, mapCgO, sc)
axis image; axis off; colormap(hsv); 
xlim([-800 800])
ylim([-800 800])
print('-depsc',[path,base,'-prop-sq']);

clf
theAlpha = -100;
ligne = find(alphaCoadd <= theAlpha, 1, 'last' );
ligne2 = find(alpha <= theAlpha, 1, 'last' );
plot(beta, mapCgO(ligne2,:) + abs(min(coaddPSW(:)/G)),'r')
hold on
plot(betaCoadd, coaddPSW(ligne,:)/G + abs(min(coaddPSW(:)/G)))
xlim([-1000 600])
print('-depsc',[path,base,'-slice']);

% %%
% figure(3)
% clf
% indexA = find(-200 <= alpha & alpha <= 100);
% indexB = find(-200 <= beta & beta <= 100);
% 
% GX = diff(mapCgO,1);
% GY = diff(mapCgO,2);
% 
% 
% subplot(221)
% imagesc(GX(indexB,indexA)); colormap(gray); axis image
% subplot(222)
% imagesc(GY(indexB,indexA)); colormap(gray); axis image
% subplot(223)
% imagesc(mapCgO(indexB,indexA)); colormap(gray); axis image
