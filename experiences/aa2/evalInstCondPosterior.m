
allvalue = linspace(3.092e4,3.112e4,100);
crit = zeros(size(allvalue));

for iIndex = 1:numel(allvalue)
  riParams{5} = allvalue(iIndex);
  [HrondProp HdirectProp] = computeRI(riParams{:});
  repro = directInvariantF(ufft2(sky), HrondProp, index250, Nalpha, Nbeta, Norder, Nbolo250, Nspeed, N_scan_total, unique_speed, the_speeds);

  dataAdeq = 0;
  for iscan = 1:N_scan_total
    dataAdeq = dataAdeq + sum(sum((data250{iscan} - repro{iscan}).^2));
  end

  priorAdeq = unifpdf(allvalue(iIndex),3.092e4,3.112e4);

  crit(iIndex) = dataAdeq;
end

figure(10000)
plot(allvalue,crit,'.')
