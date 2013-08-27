function extractAnswerkey

filename = 'myinput.m';
[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, ...
	 HSQCDATA] = loaddata(filename);
load order.m

[numPeaks,dummy] = size(HSQCDATA);


fout = fopen('myAnswerKey.m','w');
newAnswerKey = [transpose(1:numPeaks), order];
fprintf(fout, '%d %d\n', newAnswerKey');
fclose(fout);

