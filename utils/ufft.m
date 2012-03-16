function [Y]=myfft(X)
% MYFFT - The normalized Fourier transform 
%  
% * The null frequency is at 0
% * The value of the null frequency is equal to the mean*sqrt(prod(size(X)))
    
    Nx = length(X);
    
    Y = fft(X, Nx) / sqrt(Nx) ;

