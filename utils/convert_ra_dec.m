function [ xi_vect, eta_vect ] = convert_ra_dec( p,raz,decz )
%CONVERT_RA_DEC Summary of this function goes here
%   Detailed explanation goes here

%  - - - - - - - - -
% ;;   s l a D s 2 t p
% ;;  - - - - - - - - -
% ;;
% ;;  Projection of spherical coordinates onto tangent plane
% ;;  ('gnomonic' projection - 'standard coordinates').
% ;;
% ;;  (double precision)
% ;;
% ;;  Given:
% ;;     ra,dec      double   spherical coordinates of point to be projected
% ;;     raz,decz    double   spherical coordinates of tangent point
% ;;
% ;;  Returned:
% ;;     *xi,*eta    double   rectangular coordinates on tangent plane
% ;;     *j          int      status:   0 = OK, star on tangent plane
% ;;                                    1 = error, star too far from axis
% ;;                                    2 = error, antistar on tangent plane
% ;;                                    3 = error, antistar too far from axis
% ;;
% ;;  Last revision:   18 July 1996
% ;;
% ;;  Copyright P.T.Wallace.  All rights reserved.
% ;;

ra_vect=p(1,:);
dec_vect=p(2,:);

for i=1:size(ra_vect,2)
    ra=ra_vect(i);
    dec=dec_vect(i);
    
    %    disp('ra')
    %    size(ra)
    %
    %    disp('dec')
    %    size(dec)
    
    %disp('decz')
    %decz=median(dec)
    
    %disp('raz')
    %raz=median(ra)
    
    %disp('sdecz')
    sdecz = sin ( decz );
    
    %disp('sdec')
    sdec = sin ( dec );
    %size(sdec)
    
    %disp('cdecz')
    cdecz = cos ( decz );
    
    %disp('cdec')
    cdec = cos ( dec );
    %size(cdec)
    
    %disp('radif')
    radif = ra - raz;%.*ones(size(ra,1),size(ra,2));
    %size(radif)
    
    %disp('sradif')
    sradif = sin ( radif );
    %size(sradif)
    
    %disp('cradif')
    cradif = cos ( radif );
    %size(cradif)
    
    %  disp('dd etc ...')
    %  dd=sdec .* sdecz;
    %  cc=(cdec .* cdecz );
    %  size(dd)
    %  size(cc)
    %  size(cradif)
    
    denom = sdec * sdecz + (cdec * cdecz )* cradif;
    
    %size(denom)
    
    %;; Compute tangent plane coordinates (even in dubious cases)
    xi = (cdec * sradif) / denom;
    eta = ( sdec .* cdecz - (cdec * sdecz * cradif )) / denom;
    
    
    
    xi_vect(i)=xi;
    eta_vect(i)=eta;
end
%z(1,:)=xi;
%z(2,:)=eta;

% #define TINY 1e-6
% {
%    double sdecz, sdec, cdecz, cdec, radif, sradif, cradif, denom;
%
%
% /* Trig functions */
%    sdecz = sin ( decz );
%    sdec = sin ( dec );
%    cdecz = cos ( decz );
%    cdec = cos ( dec );
%    radif = ra - raz;
%    sradif = sin ( radif );
%    cradif = cos ( radif );
%
% /* Reciprocal of star vector length to tangent plane */
%    denom = sdec * sdecz + cdec * cdecz * cradif;
%
% /* Handle vectors too far from axis */
%    if ( denom > TINY ) {
%       *j = 0;
%    } else if ( denom >= 0.0 ) {
%       *j = 1;
%       denom = TINY;
%    } else if ( denom > -TINY ) {
%       *j = 2;
%       denom = -TINY;
%    } else {
%       *j = 3;
%    }
%
%    /* Compute tangent plane coordinates (even in dubious cases) */
%    *xi = cdec * sradif / denom;
%    *eta = ( sdec * cdecz - cdec * sdecz * cradif ) / denom;
% }

end