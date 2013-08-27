chRDCs = load ('C-H_medium1.m');
chRDCs(:,2) = chRDCs(:,2)*.491;
fid = fopen('C-H_medium1_reweighted.m');
for i = 1:size(chRDCs,1)
  fprintf(fid, '%f %f\n', chRDCs(i,1),chRDCs(i,2));
end
fclose(fid);
