function M = extractOrderFile

[R T RDC1 RDC2 X Y Z H_CS NH_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
    =textread('myinput.m','%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

fid = fopen('order.m','w');
for i = 1:length(R)
  fprintf(fid, '%d\n', R(i));
end

fclose(fid);