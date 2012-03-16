function x = unitifft2(image, varargin);
  %% UNITIFFT2 - The unitary inverse Fourier transform
  %%
  %% The value of the null frequency is equal to the
  %% mean*sqrt(prod(size(X)))
    
    Nx = size(image,1);  Ny = size(image,2);
    
    x = sqrt(Nx*Ny) * ifft2(image, varargin{:});
  
end
