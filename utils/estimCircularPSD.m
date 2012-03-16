function [meanPsd varargout] = estimCircularPSD(spectrum,varargin)
% ESTIMCIRCULARPSD - Compute the circular mean of the spectrum
%
% meanPsd = estimCircularPSD(SPECTRUM)
%
% return a circular mean of spectrum. The function compute the exact
% frequency axis for each elemente of spectrum and return the mean of
% frequency at the same position. The mean is sorted.
%
% Without output argument the PSD is plotted.
% 
% SPECTRUM must be 2D spectrum with the (0,0) at the upper left corner of
%            the matrix (matlab convention)
%
% meanPsd = estimCircularPSD(spectrum, DELTAF)
%
% return a circular mean of spectrum. The function compute the DELTAF
% rounded value of the frequency axis for each elemente of spectrum and
% return the mean of frequency at the same position. The mean is sorted.
%
% [m r] = estimCircularPSD(...)
% 
% return the PSD and the corresponding radius r 
%
% [m r c] = estimCircularPSD(...)
% 
% return also the number of time the frequency is in the fourrier plane (the
% mean coefficient)
% 
% $Id: estimCircularPSD.m,v 1.2 2009/05/27 11:27:29 orieux Test $
% Time-stamp: < >
    
    spectrum = abs(fftshift(spectrum)).^2;
    
% Exact frequency
    [xsize ysize] = size(spectrum);
    
    freqx = linspace(-0.5,0.5,xsize);
    freqy = linspace(-0.5,0.5,ysize);
    [FREQX,FREQY] = ndgrid(freqx,freqy);
    
% Radial frequency
    radius2D = sqrt(FREQX.^2 + FREQY.^2);
    if nargin > 1
        deltaf = varargin{1};
        radius2D = round((sqrt(FREQX.^2 + FREQY.^2))./deltaf)*deltaf;
    end
    [radius I J] = unique(radius2D(:));
    
% Donc radius2D = radius(J). Donc J contient le rayon de chaque pixel de
% psd. Je me construit une matrice pleine de zero. En ligne on a le rayon et
% donc le nombre de rayon et en colonne le nombre de pixel de psd.
    
% Avec sparse parce que sinon c'est énorme.
    formean = full(sum(sparse(J, [1:numel(spectrum)], spectrum(:), numel(radius), ...
                              numel(spectrum)), 2));
    coefs = full(sum(sparse(J, [1:numel(spectrum)], 1, numel(radius), ...
                            numel(spectrum)),2));
    
    themean = formean./coefs;
    
    if nargout == 0
        figure
        semilogy(radius,themean);
        xlabel('radius')
        return
    end
    
    meanPsd = themean;
    
    if nargout > 1
        varargout(1) = {radius};
    end
    if nargout > 2
        varargout(2) = {coefs};
    end
    
% With full matrice if you want
% formean = zeros(numel(radius), numel(spectrum));
% index = sub2ind(J,[1:numel(spectrum)]);
% formean(index) = spectrum;
% formean = sum(formean,2);
% coefs = zeros(numel(radius) numel(spectrum));
% coefs(index) = 1;
% coefs = sum(coefs,2);
% mean = formean./coefs;

% With loop
% circularMean = zeros(size(radius));
% for i = 1:length(radius)
%     pos = find(radius2D == radius(i));
%     circularMean(i) = sum(spectrum(pos)) / numel(pos);
% end

% $Log: estimCircularPSD.m,v $
% Revision 1.2  2009/05/27 11:27:29  orieux
% Use a new figure and add xlabel.
%
% Revision 1.1  2009/05/22 18:25:35  orieux
% Initial revision
%
