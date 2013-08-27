%fid = fopen('MBP_100PercentCorrect.txt','w');
%numPeaks = 335;
%numResidues = 348;
fid = fopen('1LYZ_100PercentCorrect.txt','w');
numPeaks = 126;
numResidues = 126;
for i = 1:numPeaks
  for j = 1:i-1
    fprintf(fid, ' 0 ');
  end
  fprintf(fid, ' 1 ');
  for j = i+1:numResidues
    fprintf(fid, ' 0 ');
  end
  fprintf(fid, ' %d %d %d\n', i,i,1);
end
fclose(fid);