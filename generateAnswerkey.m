name = 'InputFiles/myinput.m';

[vectors,types,resnums] = loaddata(name);

outFilename5            = 'answerkey.m';

fprintf(1, 'check out %s\n',outFilename5);

fid                     = fopen(outFilename5,'w');
    
for i = 1:length(resnums)
  fprintf(fid, '%d %d\n', i, resnums(i));
end

fclose(fid);