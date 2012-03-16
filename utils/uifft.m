function [X]=myifft(Y);
% MYIFFT - The normalized inverse Fourier transform 
%  
% * The value of the null frequency is equal to the mean*sqrt(prod(size(X)))
    
    Nx = length(Y);
    
    X = sqrt(Nx) * ifft(Y, Nx);
  
