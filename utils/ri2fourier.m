function tf = ri2fourier(ri, Nalpha, Nbeta)
  %% RI2FOURIER - compute the transfert function of the RI
  %%
  %% This function make the necessary correct zero-padding, zero
  %% convention, correct fft2 etc... to compute the transfert function
  %% of a RI.

  [Nsa Nsb] = size(ri);
        
  %% Zero padding
  ripadded = zeros(Nalpha, Nbeta);
        
  %% Place the RI before shift
  ripadded(1:Nsa,1:Nsb) = ri;
  %% Zero convention of the fft to avoid the phase problem. Work with
  %% odd and even size.
  ripadded = circshift(ripadded,[-floor(Nsa/2) -floor(Nsb/2)]);
        
  %% The fourier transforme for the impultionnal response. To use with
  %% ufft2 for signal. 
  tf = fft2(ripadded, Nalpha, Nbeta);
        
end
