dbstop if error;
dbstop if warning;
  
[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, HSQCDATA] ...
    = loaddata('myinput.m');

filename = 'parsedResonances.txt.1AAR';
fid = fopen(filename, 'w');
fprintf(1, 'check out %s\n', filename);
for i = 1:size(HSQCDATA,1)
  fprintf(fid, '%d %s %f %f\n', RESNUMS(i), TYPES{i}, HSQCDATA(i,2), HSQCDATA(i,3));
end

for i = size(HSQCDATA,1)+1:length(RESNUMS)
  fprintf(fid, '%d %s -999 -999\n', RESNUMS(i), TYPES{i});
end
fclose(fid);
