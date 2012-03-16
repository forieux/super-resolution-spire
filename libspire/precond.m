function object = precond(object, hypers, bands, Hrond250, ...
                           Hrond360, Hrond520,  coefs250, coefs360, ...
                           coefs520, regOps, objectMean, Nalpha, ...
                           Nbeta, Norder)

  compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);
  
  if compute250
    hes = 2*(hypers(1,1,1)*mean(abs(Hrond250(:,:,1,:)).^2,4) + hypers(2,1,1)*regOps(:,:,1,1));
    object(:,:,1) = object(:,:,1)./hes;
  end

  if compute360
    hes = 2*(hypers(1,2,1)*mean(abs(Hrond360(:,:,1,:)).^2,4) + hypers(2,2,1)*regOps(:,:,1,1));
    object(:,:,2) = object(:,:,2)./hes;
  end

  if compute250
    hes = 2*(hypers(1,3,1)*mean(abs(Hrond520(:,:,1,:)).^2,4) + hypers(2,3,1)*regOps(:,:,1,1));
    object(:,:,3) = object(:,:,3)./hes;
  end
  
end
