function [Hrond Hdirect] = computeRITest(lstep, bindex, lband, cl, ...
    sigmac, salpha, sbeta, astep, ...
    bstep, Nalpha, Nbeta, Nspeed, ...
    Norder, uspeed, opte, taue, gain)
%% COMPUTERI - Compute the transfert function for each band

%% (1) : Only one lambda in each band.
%% lambda = [mean(band_250') mean(band_360') mean(band_520')];

%% (2) : asymetrique

%% h(lambda)
%%  ^                          250
%%  |       *********************************************
%%  |       ^          ^          ^          ^          *
%%  |       |          |          |          |          *
%%  |       <--------->|          |          |          *
%%  |       |          |          |          |          *
%% ---------------------------------------------------------------> lambda
%%  |

%% lambda = linspace(min(lband), max(lband), lstep+1);
%% delta_lambda = lambda(2) - lambda(1);

%% (3) : symetric

%% With this option the extrem bound of the band is not actually not
%% integrated, but the integrand is symetric in respect to the central
%% wavelength

%% h(lambda)
%%  ^                          250
%%  |       *****************************************
%%  |       *        ^         ^          ^         *
%%  |       *        |         |          |         *
%%  |       *  <----------->   |          |         *
%%  |       *        |         |          |         *
%% --------------------------------------------------------------> lambda
%%  |

lambda = linspace(min(lband), max(lband), lstep+2);
lambda = lambda(2:end-1);

%% RI params
sigma = [sigmac salpha sbeta];

%% The two support must contain zero. This is the case here.
bound = round(8*sigma(1)*max(lband)./astep)*astep;
SupAlpha = -bound:astep:bound;

bound = round(8*sigma(1)*max(lband)./bstep)*bstep;
SupBeta = -bound:bstep:bound;

Nsalpha = length(SupAlpha);
Nsbeta = length(SupBeta);

Hrond = zeros(Nalpha, Nbeta, Norder, Nspeed);
Hdirect = zeros(Nsalpha, Nsbeta, Norder, Nspeed);

% uspeed=[57.475976556075466 -57.475976556075466  12.615323678396063 -12.615323678396063; ...
%     17.907282970233844 -17.907282970233844 -58.861853143721746  58.861853143721746];

for ispeed = 1:Nspeed
    
    [Hrond(:,:,:,ispeed) Hdirect(:,:,:,ispeed)] = ...
        transfertFunction(bindex, lambda, cl, sigma, ...
        uspeed(:,ispeed), Nalpha, Nbeta, Norder, ...
        SupAlpha, SupBeta, opte, taue, gain);
    
end

end

