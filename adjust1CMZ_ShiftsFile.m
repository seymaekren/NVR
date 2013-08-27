load MySHIFTS
resId = MySHIFTS(:,1);
resId = resId - 78;
fid = fopen('MySHIFTS.adjusted', 'w');
for i = 1:length(resId)
  fprintf(fid, '%d %f %f\n', resId(i), MySHIFTS(i,2), MySHIFTS(i,3));
end
fclose(fid);