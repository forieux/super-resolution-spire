function psd = psdCommonSig(frequency, frequencyCutE, declivity)
    
    psd = 1./(1 + (frequency/frequencyCutE).^declivity);
    
    psd = (psd + fliplr(psd))/1;

end
