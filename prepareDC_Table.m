%load chRdc.txt
chRdc = load ('chRdc_renumbered.txt');

%[R T RDC1 RDC2 X Y Z H_CS NH_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
%    =textread('myinput.m.hSRI.temporary','%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

[R T RDC1 RDC2 X Y Z H_CS NH_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
    =textread('myinput.m.temporary.ff2','%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

%fout = fopen('hsri_dc_newformat.tab','w');
fout = fopen('ff2_dc.tab','w');

for i = 1:size(chRdc,1)
  resId = chRdc(i,1);
  lineIndex = find(R == resId);
%  fprintf(fout, '%5d %6s %6s %5d %6s %6s %9.3f %9.3f %.2f\n', resId, ...
%  	  'deneme', 'HA', resId, 'deneme', 'CA', chRdc(i,2),1.0, ...
%  	  0.5);
  fprintf(fout, '%5d %6s %6s %5d %6s %6s %9.3f %6.3f %.2f\n', resId, ...
  	  T{lineIndex}, 'HA', resId, T{lineIndex}, 'CA', chRdc(i,2),1.0, ...
  	  0.5);
end

fprintf(fout, '\n');
%load nhRdc.txt
nhRdc = load ('nhRdc_renumbered.txt');
for i = 1:size(nhRdc,1)
  resId = nhRdc(i,1);
  lineIndex = find(R == resId);
  fprintf(fout, '%5d %6s %6s %5d %6s %6s %9.3f %6.3f %.2f\n', resId, ...
	  T{lineIndex}, 'N', resId, T{lineIndex}, 'HN', nhRdc(i,2),1.0, ...
	  1.0);
end
fclose(fout);