function arcsec = rad2arcsec(rad)
% ARCSEC2RAD - conversion de radian en arcsecondes
% 
% arcsec = rad2arcsec(rad)
% 
% rad : vecteur contenant les valeurs � convertir
% 
% arcsec : r�sultat de la conversion en arcseconde
% 
%            180*3600
%   arcsec = -------- rad;
%               pi
    
arcsec = 180/pi*3600.*rad;

