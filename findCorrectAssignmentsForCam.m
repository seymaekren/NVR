function findCorrectAssignmentsForCam


HSQC_Filename  = 'Peak-CAM13C_Dieckmann-NHSQC-freq-Refined.dat';
BMRB_Filename  = 'bmrbCAM13_3.txt.parsed';
outFilename    = 'order.m.CAM';

H_EPS = 1; N_EPS = 1;

[dummyString HSQC_N_ppm HSQC_HN_ppm intensity] = ...
    textread(HSQC_Filename, '%s %f %f %f');

HSQC_N_ppm  = round(HSQC_N_ppm*10);
HSQC_HN_ppm = round(HSQC_HN_ppm*100);

[aaNumber aaName BMRB_N_ppm BMRB_HN_ppm] = textread(BMRB_Filename, '%d %s %f %f');


BMRB_N_ppm = BMRB_N_ppm * 10;
BMRB_HN_ppm = BMRB_HN_ppm*100;

fid         = fopen(outFilename, 'w');
fprintf(1, 'check out %s\n', outFilename);


numHSQC_Peaks = length(HSQC_N_ppm);

numUnassigned = 0;
for i = 1:numHSQC_Peaks
  bmrbIndex = findCorrespondingBMRB_Peak(HSQC_N_ppm(i), HSQC_HN_ppm(i), ...
					BMRB_N_ppm, ...
					BMRB_HN_ppm, H_EPS, N_EPS);
  if (bmrbIndex ~= -1)
    resIndex  = aaNumber(bmrbIndex);
  else
    resIndex  = -1;
    numUnassigned = numUnassigned + 1;
  end
  fprintf(fid, '%d %d\n', i, resIndex);
end

fprintf(1, '%d peaks out of %d have been unassigned.\n',numUnassigned, ...
	numHSQC_Peaks);


%this function          returns the bmrb peak index corresponding to the
%HSQC Peak.
function bmrbPeakIndex = findCorrespondingBMRB_Peak(hsqcN_ppm, hsqc_HN_ppm,...
						  BMRB_N_ppm, BMRB_HN_ppm,...
						  H_EPS, ...
						  N_EPS)


closeBMRB_Peaks = [];
numBMRB_Peaks = length(BMRB_N_ppm);

for i = 1:numBMRB_Peaks
  N_CS_Diff = abs(hsqcN_ppm-BMRB_N_ppm(i));
  H_CS_Diff = abs(hsqc_HN_ppm-BMRB_HN_ppm(i));
% if ((H_CS_Diff < H_EPS) & (N_CS_Diff < N_EPS))
  if ((H_CS_Diff == 0) & (N_CS_Diff == 0))
     closeBMRB_Peaks = [closeBMRB_Peaks i];
  end
end

if (length(closeBMRB_Peaks) == 1)
  bmrbPeakIndex = closeBMRB_Peaks;
elseif (length(closeBMRB_Peaks) > 1)
  fprintf(1, 'this HSQC peak has more than one BMRB peak close to it.\n');
  closeBMRB_Peaks
  hsqcN_ppm
  hsqc_HN_ppm
  BMRB_N_ppm(closeBMRB_Peaks)
  BMRB_HN_ppm(closeBMRB_Peaks)
  keyboard
  bmrbPeakIndex = -1;
elseif (length(closeBMRB_Peaks) == 0)
  bmrbPeakIndex = -1;
end
