function [meanPsd varargout] = estimCircularPSDradon(im,varargin)
% ESTIMCIRCULARPSDRADON - Compute the circular mean of the spectrum with
% radon transform
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
% $Id: estimCircularPSDRadon.m,v 1.1 2009/07/22 17:49:58 orieux Exp $
% Time-stamp: < >
    
%% PSD for each angle
    spectrum = abs(fft(radon(im))).^2;
    %% Radial symetry and mean on angle
    spectrum = mean(spectrum' + flipud(spectrum)');
    spectrum = spectrum(1:floor(end/2));
    
    if nargout == 0
        figure
        semilogy(radius,spectrum);
        xlabel('radius')
        return
    end
    
    meanPsd = spectrum;
    
    if nargout > 1
        varargout(1) = {radius};
    end
    
% $Log: estimCircularPSDRadon.m,v $
% Revision 1.1  2009/07/22 17:49:58  orieux
% Initial revision
%
