function rearrangeCamData


%this function will put the correct chemical shifts for each line
%in myinput.m for CAM. If the correct chemical shifts for an amino
%acid is unknown, -999 will be put there. This way we will have
%less peaks to assign but we will know the assignments of those peaks.

outfilename = 'myinput.m.CAM.rearranged';
outfilename2 = 'myinput.m.CAM.chemicalShiftsShuffled';
fprintf(1, 'check out %s and %s\n', outfilename,outfilename2);

fid = fopen(outfilename2, 'w');

[R T RDC1 RDC2 X Y Z H_CS N_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
    =textread('myinput.m','%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

ORDER = load ('order.m.CAM');

h_ppm = ones(max(R),1)*-999.0;
n_ppm = h_ppm;

for i = 1:length(ORDER)
  peakIndex    = ORDER(i,1);
  residueIndex = ORDER(i,2);
  if (residueIndex ~= -1)
    assert (peakIndex == i);
    h_ppm(residueIndex) = H_CS(i);
    n_ppm(residueIndex) = N_CS(i);
  end
end

for i = 1:length(R)
  aaIndex   = R(i);
  if (aaIndex == -999)
    continue;
  end
  myH_CS    = h_ppm(aaIndex);
  myN_CS    = n_ppm(aaIndex);
  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%d\t%f\t%f\t%f\n', R(i), T{i}, RDC1(i), RDC2(i), X(i),Y(i),Z(i),myH_CS, myN_CS,SS{i},HB{i},MOLMOL_HB(i),HX(i),HY(i),HZ(i));
end

fclose(fid);

[R T RDC1 RDC2 X Y Z H_CS N_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
    =textread(outfilename2,'%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

fid = fopen(outfilename, 'w');

printLast = [];
for i = 1:length(R)
  if (H_CS(i) == -999)
    printLast = [printLast i];
  else
    fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%d\t%f\t%f\t%f\n', R(i), T{i}, RDC1(i), RDC2(i), X(i),Y(i),Z(i),H_CS(i), N_CS(i),SS{i},HB{i},MOLMOL_HB(i),HX(i),HY(i),HZ(i));
  end
end

for j = 1:length(printLast)
  i = printLast(j);
  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%d\t%f\t%f\t%f\n', R(i), T{i}, RDC1(i), RDC2(i), X(i),Y(i),Z(i),H_CS(i), N_CS(i),SS{i},HB{i},MOLMOL_HB(i),HX(i),HY(i),HZ(i));
end
fclose(fid);