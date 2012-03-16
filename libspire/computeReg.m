function [diffAlpha diffBeta circDalpha circDbeta] = ...
      computeReg(salpha, sbeta, skya, skyb, N, Nalpha, Nbeta)
  %% COMPUTEREG - Compute the regularisation operator
  
  diffAlpha = diffOpAlpha(salpha, sbeta, skya, skyb, N);
  diffBeta = diffOpBeta(salpha, sbeta, skya, skyb, N);

  circDalpha = real(ri2fourier(diffAlpha, Nalpha, Nbeta));
  circDbeta  = real(ri2fourier(diffBeta, Nalpha, Nbeta));

end

