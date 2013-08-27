function printTen(S, filename)

fid = fopen(filename, 'w');
fprintf(1, 'check out %s\n', filename);
assert (size(S,1) == 3);
assert (size(S,2) == 3);
for i = 1:3
  for j = 1:3
    fprintf(fid, '%d ', S(i,j));
  end
  fprintf(fid, '\n');
end
fclose(fid);