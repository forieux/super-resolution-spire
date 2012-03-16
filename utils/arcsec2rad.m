function rad = arcsec2rad(arcsec)
% ARCSEC2RAD - conversion de arc-seconde en radian
%   
% rad = arcsec2rad(arcsec)
% 
% arcsec : vecteur de valeur en arc-secondes
% 
% rad : résultat de la conversion en radian
%
%             pi     
%    rad = -------- arcseconde
%          180*3600
  
rad = (arcsec./3600).*pi./180;

