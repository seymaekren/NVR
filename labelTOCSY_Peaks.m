function labelTOCSY_Peaks

TOCSY_Filename  = 'Peak-CAM13C_Dieckmann-TOCSYHSQC-freq-Refined.dat';
HSQC_Filename   = 'Peak-CAM13C_Dieckmann-NHSQC-freq-Refined.dat';
outFilename     = 'labeledTOCSY_Peaks.txt';

H_EPS = 0.03; N_EPS = 0.3;

[dummyString TOCSY_H_sidechain_ppm TOCSY_N1_ppm TOCSY_H1_ppm intensity] = ...
    textread(TOCSY_Filename, '%s %f %f %f %f');

[dummyString HSQC_N_ppm HSQC_HN_ppm intensity] = ...
    textread(HSQC_Filename, '%s %f %f %f');

numTOCSY_Peaks  = length(TOCSY_H1_ppm);

fid         = fopen(outFilename, 'w');
fprintf(1, 'check out %s\n', outFilename);


for tocsyPeakIndex = 1:numTOCSY_Peaks
  [closeHSQC_PeakIndices] = findCloseHSQC_Peak(TOCSY_N1_ppm(tocsyPeakIndex), ...
				      TOCSY_H1_ppm(tocsyPeakIndex), ...
				      HSQC_N_ppm, HSQC_HN_ppm, H_EPS, ...
					       N_EPS);
  if (length(closeHSQC_PeakIndices) == 1)
    if (TOCSY_H_sidechain_ppm(tocsyPeakIndex) > 0)
      %not printing negative chemical shifts, don't know what they mean?
      fprintf(fid, '%f %f 100 %d\n', TOCSY_H_sidechain_ppm(tocsyPeakIndex), ...
	      TOCSY_H_sidechain_ppm(tocsyPeakIndex), ...
	      closeHSQC_PeakIndices(1));
    end
  end
end

fprintf(1, 'filtered out negative chemical shifts.\n');

function [closeHSQC_PeakIndices] = findCloseHSQC_Peak(TOCSY_N_ppm, ...
						  TOCSY_H1_ppm, ...
						  HSQC_N_ppm, HSQC_HN_ppm, ...
						  H_EPS, N_EPS);

numHSQCPeaks = size(HSQC_N_ppm,1);
closeHSQC_PeakIndices = [];
for peakIndex = 1:numHSQCPeaks
  N_CS_DIFF  = abs(TOCSY_N_ppm-HSQC_N_ppm(peakIndex));
  H_CS_DIFF  = abs(TOCSY_H1_ppm-HSQC_HN_ppm(peakIndex));
  
  if ((N_CS_DIFF < N_EPS) & (H_CS_DIFF < H_EPS))
    closeHSQC_PeakIndices = [closeHSQC_PeakIndices peakIndex];
  end
end


