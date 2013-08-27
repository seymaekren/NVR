function hnClosePeaks = findCloseToHN(h_ppm, n_ppm, HSQCDATA, HN_EPS, ...
				      N_EPS)

%finds the HSQC peak indices that are close to the input chemical
%shifts h_ppm and n_ppm.


%HN_EPS = 0.03; N_EPS = 0.3; %these thresholds are also
			    %set in findCloseToHN and findCloseToH and findUnambiguousNOEs.m.
numPeaks = size(HSQCDATA,1);
hnClosePeaks = [];

%fprintf(1, 'the input chemical shifts are: %f %f\n', h_ppm, n_ppm);

for peakIndex = 1:numPeaks
  H_CS      = HSQCDATA(peakIndex,2);
  H_CS_Diff = abs(H_CS - h_ppm);
  N_CS      = HSQCDATA(peakIndex,3);
  N_CS_Diff = abs(N_CS - n_ppm);
  if (H_CS_Diff < HN_EPS) & (N_CS_Diff < N_EPS)
    hnClosePeaks = [hnClosePeaks peakIndex];
  end
end
