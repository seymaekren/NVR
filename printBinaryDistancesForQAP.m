function printBinaryDistancesForQAP(COLIN, NTH, ALLDISTS);


DISTANCE_THRESHOLD = 6;
fid = fopen('BinaryDistancesForQAP.txt','w');
fprintf(1, 'printing to BinaryDistancesForQAP.txt\n');
for residue1Index = 1:length(COLIN)
  for residue2Index =  1:length(COLIN)
    if (ALLDISTS(COLIN(residue1Index),COLIN(residue2Index)) > DISTANCE_THRESHOLD)
      fprintf(fid, '%f ', ALLDISTS(COLIN(residue1Index), ...
				   COLIN(residue2Index)) - DISTANCE_THRESHOLD);
    else
      fprintf(fid, '0 ');
    end
  end	
  fprintf(fid, '\n');
end
fclose(fid);
