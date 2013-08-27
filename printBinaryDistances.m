function printBinaryDistances(COLIN, NTH, ALLDISTS);

fid = fopen('BinaryDistances.txt','w');
fprintf(1, 'printing to BinaryDistances.txt\n');
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
