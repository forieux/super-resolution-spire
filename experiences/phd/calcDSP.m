cielRond = zeros(size(map));
for iband=1:3
    cielRond(:,:,iband) = ufft2(map(:,:,iband));
end

cielRondMean = zeros(size(mapmean));
for iband=1:3
    cielRondMean(:,:,iband) = ufft2(mapmean(:,:,iband));
end

% mapinterp = zeros(size(mapnan));
% mapinterp_dspest = zeros(size(mapnan));
% cielRondInterp = zeros(size(map));
% for iband = 1:3
%     [I J] = ind2sub([Nalpha Nbeta], find(isnan(mapnan(:,:,1)) == 0));
%     mapinterp(:,:,iband) = detNan(mapnan(:,:,iband));
% % % Remove extrapolation for better comparison of the mean
% %     mapinterp_dspest(min(I):max(I),min(J):max(J),iband) = mapinterp(min(I):max(I),min(J):max(J),iband);
%     cielRondInterp(:,:,iband) = ufft2(mapinterp(:,:,iband));
% end

cielRondMADmap = zeros(size(madmap));
for iband = 1:3,
    cielRondMADmap(:,:,iband) = ufft2(madmap(:,:,iband));
end

%% MADmap (dirty map with 6 arcsec resolution)


%% Computation of DSP
deltaF = 0.005;
[dspmap250 radialFreq] = estimCircularPSD(ufft2(map(:,:,1)/Hrond250(1,1,1)),deltaF);
dspmap360 = estimCircularPSD(ufft2(map(:,:,2)/Hrond360(1,1,1)),deltaF);
dspmap520 = estimCircularPSD(ufft2(map(:,:,3)/Hrond520(1,1,1)),deltaF);

% dspmapinterp250 = estimCircularPSD(ufft2(mapinterp_dspest(:,:,1)./Hrond250(1,1,1)),deltaF);
% dspmapinterp360 = estimCircularPSD(ufft2(mapinterp_dspest(:,:,2)./Hrond360(1,1,1)),deltaF);
% dspmapinterp520 = estimCircularPSD(ufft2(mapinterp_dspest(:,:,3)./Hrond520(1,1,1)),deltaF);

[dspmadmap250 radialFreq_mm]= estimCircularPSD(ufft2(madmap(:,:,1)./ ...
                                                  Hrond250(1,1,1)),deltaF);
dspmadmap360 = estimCircularPSD(ufft2(madmap(:,:,2)./Hrond360(1,1,1)),deltaF);
dspmadmap520 = estimCircularPSD(ufft2(madmap(:,:,3)./Hrond520(1,1,1)),deltaF);

dspnoisefree360 = estimCircularPSD(ufft2(sky(:,:,2)).*Hrond360(:,:,1)./ ...
                                   Hrond360(1,1,1),deltaF);

dspsky250 = estimCircularPSD(ufft2(sky(:,:,1)),deltaF);
dspsky360 = estimCircularPSD(ufft2(sky(:,:,2)),deltaF);
dspsky520 = estimCircularPSD(ufft2(sky(:,:,3)),deltaF);

% To compare dirty, mine and madmap, when need to fixe the absolute freq,
% not use only reduced freq, since the sampling is not the same.
limitFreq = 0.03; % arcsec^-1
Nfreq = length(radialFreq);
Te = mean([alpha_step beta_step]);
maxFreq = max(radialFreq);
indFreq = round(limitFreq*Nfreq*Te/maxFreq);

Nfreq = length(radialFreq_mm);
Te = 6;
maxFreq = max(radialFreq_mm);
indFreq_mm = round(limitFreq*Nfreq*Te/maxFreq);

radialFreq_cut = radialFreq(1:indFreq);
radialFreq_mm_cut = radialFreq_mm(1:indFreq_mm);

dspmap250_cut  = dspmap250(1:indFreq);
dspmap360_cut  = dspmap360(1:indFreq);
dspmap520_cut  = dspmap520(1:indFreq);

% dspmapinterp250_cut  = dspmapinterp250(1:indFreq);
% dspmapinterp360_cut  = dspmapinterp360(1:indFreq);
% dspmapinterp520_cut  = dspmapinterp520(1:indFreq);

dspmadmap250_cut  = dspmadmap250(1:indFreq_mm);
dspmadmap360_cut  = dspmadmap360(1:indFreq_mm);
dspmadmap520_cut  = dspmadmap520(1:indFreq_mm);

dspnoisefree360_cut  = dspnoisefree360(1:indFreq_mm);

dspsky250_cut  = dspsky250(1:indFreq);
dspsky360_cut  = dspsky360(1:indFreq);
dspsky520_cut  = dspsky520(1:indFreq);
