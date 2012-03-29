function observed = calcObsPix(coefs, Nalpha, Nbeta, Nscan, size)

  observed = zeros(Nalpha, Nbeta);
  for iscan = 1:Nscan
    observed(find(coefs{iscan} ~= 0)) = 1;
  end
  % observed = imerode(imdilate(observed, strel('diamond', 1)), strel('diamond', 1));

  % observed = imerode(observed, strel('diamond', size));

end
