order= load ('order.m');
NOEs = load ('NOEs.txt.CAM2');
noes = convertResidueIndicesToPeakIndicesInNOE_File(NOEs, order);
filename = 'NOEs_peakIndices.txt.CAM2';
fprintf(1, 'check out %s\n', filename);
fid = fopen(filename, 'w');
for i = 1:size(noes,1)
  fprintf(fid, '%d %d\n',noes(i,1),noes(i,2));
end
fclose(fid);