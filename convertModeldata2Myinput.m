function convertModeldata2Myinput

[RESNUMS AATYPE X Y Z  SS, HB,HX,HY,HZ] =textread('Mode7.coords6.1EF1','%f %s  %f %f %f %s %s %f %f %f');

HSQCDATA = load('hsqcdata.m.1EF1');

outfilename = 'myinput.m.1EF1.new';
fprintf(1, 'check out %s\n', outfilename);
fid         = fopen(outfilename, 'w');

for i = 1:size(HSQCDATA,1)

  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%d\t%f\t%f\t%f\n', RESNUMS(i), AATYPE{i}, HSQCDATA(i,4), HSQCDATA(i,5), ...
	  X(i),Y(i),Z(i),HSQCDATA(i,2), HSQCDATA(i,3), ...
	  SS{i}, HB{i}, HSQCDATA(i,6), HX(i),HY(i),HZ(i));
end

fclose(fid);

fprintf(1, 'enter return to exit debug mode.\n');
keyboard
		

