[h_ppm n_ppm h_ppm2 intensity] = ...
    textread('hnNoeManualPeak.txt.NOE.withIntensities','%f %f %f %f');
dnnIndices = [];

for crpIndex = 1:length(h_ppm)
  if (h_ppm2(crpIndex) > 6)
    dnnIndices = [dnnIndices crpIndex];
  end
end

fid = fopen('dnns.txt','w');
fprintf(1, 'check out dnns.txt\n');
for i = 1:length(dnnIndices)
  fprintf(fid, '%f %f %f %f\n',h_ppm(dnnIndices(i)),n_ppm(dnnIndices(i)), h_ppm2(dnnIndices(i)),intensity(dnnIndices(i)));
end
fclose(fid);


