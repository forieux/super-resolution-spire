function [Hrond Hdirect] = transfertFunction(band, lambda, lambda0, sigma, ...
                                             speed, Nalpha, Nbeta, Norder, ...
                                             SupAlpha, SupBeta, opte, taue, gain)
%% TRANSFERTFUNCTION - Compute the transfert function for one band at all
%%
%% [Hrond Hdirect] = transfertFunction(lambda, lambda0, sigma, speed, Nalpha,
%% Nbeta, Norder, SupAlpha, SupBeta)
%%                                  
%% compute the transfert function, and impultionnal response, for a band for
%% all order for one speed. Make the lambda numerical integration.
%% 
%% INPUT PARAMETERS
%% 
%% band -- the band 1 for 250, 2 for 260 and 3 for 520. Necessary to know the
%% filter
%% 
%% lambda -- the lambda discretisation for the integration
%% 
%% lambda0 -- the central wavelength for the lambda decomposition
%% 
%% sigma -- The vector that define the dependence in lambda.
%% 
%% speed -- The vector speed
%% 
%% Nalpha, Nbeta -- the number of alpha and beta
%% 
%% Nordre -- the number of order for lambda decomposition
%% 
%% SupAlpha, SupBeta -- the support on alpha and beta to use for the IR
%% computation
    
%% See ri_fullgaussian for other description
 
    delta_lambda = lambda(2) - lambda(1);
    
    Nsalpha = length(SupAlpha);
    Nsbeta = length(SupBeta);
    
    filterValue = lambdaFilter(band, lambda);
    
    Hdirect = zeros(Nsalpha, Nsbeta, Norder);
    Hrond = zeros(Nalpha, Nbeta, Norder);
  
    %% Computation of ri for each order for each band for each speed
    for iorder = 0:Norder-1
        
        H = ri_fullgaussian(SupAlpha, SupBeta, lambda, opte, ...
                            sigma, taue, gain, speed);
        
        %% Lambda ponderation for iorder
        pond = filterValue.*delta_lambda.*(lambda - lambda0).^iorder;
        %% use reshape+repmat because repmat alone doesn't work
        pond = reshape(pond, 1, 1, length(lambda));
        pond = repmat(pond, [Nsalpha Nsbeta 1]);
        %% sum on lambda
        Hsum = sum(H.*pond,3);
        Hdirect(:,:,iorder+1) = Hsum;
        Hrond(:,:,iorder+1) = ri2fourier(Hsum, Nalpha, Nbeta);
        
        %% Uncomment next line if you want to compute the transpose.
        %Hrond(:,:,iorder+1)= ri2fourier(fliplr(flipud(Hp)));
    
    end
    
end
