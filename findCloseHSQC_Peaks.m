%this computation is for EIN imo 4/10/2010

[residueID residueName H_CS N_CS X Y Z] = textread('combinedResonancesAndProtonCoordinates.txt','%d %s %f %f %f %f %f');
%HN_ChemicalShifts = load ('HN_chemicalShifts.txt');

%HSQC_Filename  = 'Peak-CAM13C_Dieckmann-NHSQC-freq-Refined.dat';

%[dummyString N_CS H_CS intensity] = ...
%    textread(HSQC_Filename, '%s %f %f %f');

%residueID   = HN_ChemicalShifts(:,1);
%N_CS        = HN_ChemicalShifts(:,2);
%H_CS        = HN_ChemicalShifts(:,3);

deltaH   = 0.02;
deltaN   = 0.2;
filename = 'ambiguousHSQC_ResidueIndices.txt';
fid      = fopen(filename,'w');
fprintf(1, 'check out %s\n', filename);

for i = 1:length(H_CS)
  for j = i+1:length(H_CS)
    if (abs(H_CS(i)-H_CS(j)) <= deltaH) & (abs(N_CS(j)-N_CS(i)) <= deltaN)
      fprintf(1,'entry %d and %d are ambiguous HSQC peaks.\n',i,j);
      fprintf(1, 'corresponds to residues %d and %d\n', residueID(i), ...
	      residueID(j));
      fprintf(fid, '%d %d\n',residueID(i),residueID(j));
    end
  end
end

fclose(fid);