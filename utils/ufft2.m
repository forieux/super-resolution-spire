function y = unitfft2(image, varargin)
  %% UNITFFT2 - The unitary Fourier transform
  %%
  %% The value of the null frequency is equal to the
  %% mean*sqrt(prod(size(IMAGE)))
    
    Nx = size(image,1); Ny = size(image,2);
    
    y = fft2(image, varargin{:}) / sqrt(Nx*Ny);

end
