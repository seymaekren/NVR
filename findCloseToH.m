function hClosePeaks = findCloseToH(h_ppm, HSQCDATA, H_EPS, N_EPS)

%H_EPS        = 0.03; %these thresholds are also
		     %set in findCloseToHN and
		     %findCloseToH and findUnambiguousNOEs.m
		     
numPeaks     = size(HSQCDATA,1);
hClosePeaks  = [];

for peakIndex   = 1:numPeaks
  H_CS          = HSQCDATA(peakIndex,2);
  H_CS_Diff     = abs(H_CS - h_ppm);
  if (H_CS_Diff < H_EPS) 
    hClosePeaks = [hClosePeaks peakIndex];
  end
end
