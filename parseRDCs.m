dbstop if error;
dbstop if warning;
  
[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, HSQCDATA] ...
    = loaddata('myinput.m');

filename = 'nhRdc.m.1AAR.1';
fid = fopen(filename, 'w');
fprintf(1, 'check out %s\n', filename);
for i = 1:size(HSQCDATA,1)
  fprintf(fid, '%d %f  0.0000  0.0000\n', RESNUMS(i),  HSQCDATA(i,4));
end

fclose(fid);

filename = 'nhRdc.m.1AAR.2';
fid = fopen(filename, 'w');
fprintf(1, 'check out %s\n', filename);
for i = 1:size(HSQCDATA,1)
  fprintf(fid, '%d %f  0.0000  0.0000\n', RESNUMS(i),  HSQCDATA(i,5));
end

fclose(fid);
