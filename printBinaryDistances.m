function printBinaryDistances(COLIN, NTH, ALLDISTS)

filename = sprintf('OutputFiles/BinaryDistances.txt');
fid = fopen(filename,'w');
fprintf(1, 'printing to %s\n', filename);
for residue1Index = 1:length(COLIN)
  for residue2Index =  1:length(COLIN)
    if (ALLDISTS(COLIN(residue1Index),COLIN(residue2Index)) > NTH)
      fprintf(fid, '0 ');
    else
      fprintf(fid, '1 ');
    end
  end	
  fprintf(fid, '\n');
end
fclose(fid);
