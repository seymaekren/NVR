function printDistances(ALLDISTS)
fid = fopen('distances.txt', 'w');
fprintf(1, 'check out distances.txt\n');
[m,n] = size(ALLDISTS);
for i = 1:m
  for j = 1:n
    fprintf(fid, '%f ', ALLDISTS(i,j));
  end
  fprintf(fid, '\n');
end
fclose(fid);