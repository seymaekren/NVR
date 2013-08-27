filename   = '1EZA.NOE_List.HnHn';

parsedNOEs = load (filename);

outputFilename = '1EZA.NOE_List.HnHn.C_TerminusRemoved';

fid = fopen(outputFilename,'w');
fprintf(1, 'check out %s\n', outputFilename);

for i = 1:size(parsedNOEs,1)
  if ((parsedNOEs(i,1) >= 250) & (parsedNOEs(i,1) <= 259))
    continue;
  end
  if ((parsedNOEs(i,2) >= 250) & (parsedNOEs(i,2) <= 259))
    continue;
  end
  fprintf(fid, '%d %d\n',parsedNOEs(i,1),parsedNOEs(i,2));
end