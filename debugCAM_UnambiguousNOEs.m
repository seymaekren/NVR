function debugCAM_UnambiguousNOEs

NOESY_Filename = 'Peak-CAM13C_Dieckmann-NNOESY-freq-Refined.dat';
HSQC_Filename  = 'Peak-CAM13C_Dieckmann-NHSQC-freq-Refined.dat';


H_EPS = 0.03; N_EPS = 0.3;

%peak1Index = 119; peak2Index = 126;
peak1Index = 37; peak2Index = 148;
%there exists an NOE between these two HSQC peaks that is
%unambiguous.
%there probably exists two such NOEs, one between peak1 and peak2,
%and the other, btw. peak2 and peak1.


[dummyString NOE_HN_ppm NOE_N_ppm NOE_H_ppm intensity] = ...
    textread(NOESY_Filename, '%s %f %f %f %f');

[dummyString HSQC_N_ppm HSQC_HN_ppm intensity] = ...
    textread(HSQC_Filename, '%s %f %f %f');

numNOEs  = length(NOE_HN_ppm);

HSQCDATA = zeros(length(HSQC_N_ppm), 3);

for i = 1:length(HSQC_N_ppm)
  HSQCDATA(i,3) = HSQC_N_ppm(i);
  HSQCDATA(i,2) = HSQC_HN_ppm(i);
end



n1 = HSQC_N_ppm(peak1Index)
h1 = HSQC_HN_ppm(peak1Index)
n2 = HSQC_N_ppm(peak2Index)
h2 = HSQC_HN_ppm(peak2Index)

NOEs = zeros(numNOEs, 3);
for i = 1:numNOEs
  NOEs(i, 1) = NOE_HN_ppm(i);
  NOEs(i, 2) = NOE_N_ppm (i);
  NOEs(i, 3) = NOE_H_ppm (i);
end

noeIndices = [];

for noeIndex = 1:numNOEs
  if (chemicalShiftsClose(NOEs(noeIndex,1),h1, H_EPS) & chemicalShiftsClose(NOEs(noeIndex,2),n1,N_EPS) & ...
      chemicalShiftsClose(NOEs(noeIndex,3),h2, H_EPS))
    noeIndices = [noeIndices noeIndex];
  elseif (chemicalShiftsClose(NOEs(noeIndex,3),h1, H_EPS) & chemicalShiftsClose(NOEs(noeIndex,2),n1, N_EPS) & ...
	  chemicalShiftsClose(NOEs(noeIndex,1),h2, H_EPS))
    noeIndices = [noeIndices noeIndex];
  elseif (chemicalShiftsClose(NOEs(noeIndex,1),h2, H_EPS) & chemicalShiftsClose(NOEs(noeIndex,2),n2, N_EPS) & ...
      chemicalShiftsClose(NOEs(noeIndex,3),h1, H_EPS))
    noeIndices = [noeIndices noeIndex];
  elseif (chemicalShiftsClose(NOEs(noeIndex,3),h2, H_EPS) & chemicalShiftsClose(NOEs(noeIndex,2),n2, N_EPS) & ...
	  chemicalShiftsClose(NOEs(noeIndex,1),h1, H_EPS))
    noeIndices = [noeIndices noeIndex];
  end
end
fprintf(1, 'the noe indices matching the noe between %d and %d are:\n',peak1Index,peak2Index);
noeIndices
for i = 1:length(noeIndices)
  fprintf(1, 'NOE %d: %f %f %f\n', noeIndices(i),NOEs(noeIndices(i),1), NOEs(noeIndices(i),2), NOEs(noeIndices(i),3));
end

function retVal = chemicalShiftsClose(v1,v2,eps)
if (abs(v1-v2)<eps)
  retVal = 1;
else
  retVal = 0;
end
