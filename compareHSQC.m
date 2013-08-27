[extrN pnt1N ppm1 pnt2 ppm2  magn] = ...
    textread('hsqc.lis','%d %f %f %f %f %f');

[R T RDC1 RDC2 X Y Z H_CS NH_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
    =textread('myinput.m','%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

H_EPS = 0.005;
N_EPS = 0.05;
foundPeaksInMyHSQC   = zeros(length(R),1);
foundPeaksInAnsig    = zeros(length(ppm1),1);

for i = 1:length(ppm1)
  numFoundCloseHSQC_Peaks = 0;
  for j = 1:length(H_CS)
    if (abs(ppm1(i) - H_CS(j)) < H_EPS) & (abs(ppm2(i) - NH_CS(j)) < N_EPS)
      numFoundCloseHSQC_Peaks = numFoundCloseHSQC_Peaks + 1;
      if (numFoundCloseHSQC_Peaks > 1)
	break;
      end

      peakIndex    = j;
      residueIndex = R(j);
    end
  end
  if (numFoundCloseHSQC_Peaks > 1)
    fprintf(1, 'found multiple closeby residues to %d th hsqc\n',i);
  elseif (numFoundCloseHSQC_Peaks == 1)
    fprintf(1, 'found %d th hsqc from the list to be close to',i);
    fprintf(1, ' residue %d\n', residueIndex);
    foundPeaksInAnsig(i)          = 1;
    foundPeaksInMyHSQC(peakIndex) = 1;
  end
end

unassignedPeakIndicesInMyHSQC = find(foundPeaksInMyHSQC == 0);
fprintf(1, 'unassigned peak indices in my hsqc are:\n')
for (i = 1:length(unassignedPeakIndicesInMyHSQC))
  fprintf(1, '%d,',unassignedPeakIndicesInMyHSQC(i));
end
fprintf(1, '\n');

fprintf(1, 'extra lines that could not be assigned in ANSIG file:');
extraHSQC_Peaks_In_Ansig_File = find(foundPeaksInAnsig == 0);
for i = 1:length(extraHSQC_Peaks_In_Ansig_File)
  fprintf(1, '%d,',extraHSQC_Peaks_In_Ansig_File(i));
end
fprintf(1, '\n');