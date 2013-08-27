fid = fopen('optimumScoringMatrix.txt','w');
for i = 1:70
  for j = 1:72
    if (i == j)
      fprintf(fid, '%f ', 0);
    else
      fprintf(fid, '%f ', 1E+9);
    end
  end
  fprintf(fid, '\n');
end