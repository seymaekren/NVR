function parseXpkFile
NOESY_Filename = 'ff2_n15_sansHeader_modified.xpk';
%NOESY_Filename = 'deneme.txt';
[noeIndex dummyString NOE_HN_ppm dummyValue dummyValue dummyString ...
 dummyValue dummyString dummyString NOE_H_ppm dummyValue ...
 dummyValue dummyString dummyValue dummyString dummyString ...
 NOE_N_ppm dummyValue dummyValue dummyString dummyValue ...
 dummyString dummyValue dummyValue dummyValue ...
 dummyString dummyValue] =     textread(...
     NOESY_Filename, ...
     '%d %s %f %f %f %s %f %s %s %f %f %f %s %f %s %s %f %f %f %s %f %s %f %f %f %s %f');

numNOEs  = length(NOE_HN_ppm);

NOEs = zeros(numNOEs, 3);
for i = 1:numNOEs
  NOEs(i, 1) = NOE_HN_ppm(i);
  NOEs(i, 2) = NOE_N_ppm (i);
  NOEs(i, 3) = NOE_H_ppm (i);
end


fid = fopen('ff2_n15_extractedNOEs.txt', 'w');
fprintf    (1, 'check out ff2_n15_extractedNOEs.txt\n');

for i = 1:numNOEs
  fprintf(fid, '%f %f %f\n', NOEs(i,1), NOEs(i,2), NOEs(i,3));
end