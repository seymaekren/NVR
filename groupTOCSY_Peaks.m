function groupTOCSY_Peaks

TOCSY_Filename  = 'Peak-CAM13C_Dieckmann-TOCSYHSQC-freq-Refined.dat';
outFilename    = 'groupedTOCSY_Peaks.txt';

H_EPS = 1E-6; N_EPS = 1E-6;

[dummyString H_sidechain_ppm N1_ppm H1_ppm intensity] = ...
    textread(TOCSY_Filename, '%s %f %f %f %f');

numTOCSY_Peaks  = length(H1_ppm);

fid         = fopen(outFilename, 'w');
fprintf(1, 'check out %s\n', outFilename);


groupedPeaks = zeros(numTOCSY_Peaks,1);

for peak1Index = 1:numTOCSY_Peaks

  closePeaks = [peak1Index];
  
  for peak2Index = peak1Index+1:numTOCSY_Peaks
    
    if (peaksClose(peak1Index,peak2Index,N1_ppm, ...
		   H1_ppm, H_EPS, N_EPS))
  
      closePeaks = [closePeaks peak2Index];
    end
  end

  if (length(closePeaks)>1)
  
    printClosePeaks(fid, closePeaks, H_sidechain_ppm, N1_ppm, H1_ppm);

    for i = 1:length(closePeaks)
      groupedPeaks(closePeaks(i)) = 1;
    end
  end
    
end

printUngroupedPeaks(fid,groupedPeaks,H_sidechain_ppm, N1_ppm, ...
			 H1_ppm);


function retval = peaksClose(peak1Index,peak2Index,N1_ppm, ...
			     H1_ppm, H_EPS, N_EPS)

retval = 0;
N_CS_DIFF  = abs(N1_ppm(peak1Index)-N1_ppm(peak2Index));
H_CS_DIFF  = abs(H1_ppm(peak1Index)-H1_ppm(peak2Index));

if ((N_CS_DIFF < N_EPS) & (H_CS_DIFF < H_EPS))
  retval            = 1;
end

function printClosePeaks(fid, closePeaks, H_sidechain_ppm, N1_ppm, ...
			 H1_ppm);
for i = 1:length(closePeaks)
  fprintf(fid, '%d ',closePeaks(i));
end
fprintf(fid, '\n');

for i = 1:length(closePeaks)
  fprintf(fid, '%f %f %f\n ',H_sidechain_ppm(closePeaks(i)), N1_ppm(closePeaks(i)), ...
	  H1_ppm(closePeaks(i)));
end
fprintf(fid, '\n\n');

function printUngroupedPeaks(fid,groupedPeaks,H_sidechain_ppm, N1_ppm, ...
			     H1_ppm);

ungroupedPeakIndices = find(groupedPeaks == 0);
for i = 1:length(ungroupedPeakIndices)
  fprintf(fid, '%d\n', ungroupedPeakIndices(i));
  fprintf(fid, '%f %f %f\n',H_sidechain_ppm(ungroupedPeakIndices(i)), N1_ppm(ungroupedPeakIndices(i)), ...
	  H1_ppm(ungroupedPeakIndices(i)));
end