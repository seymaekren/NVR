function M = stripTOCSY_File

[TH1,TH2,TN,TRN,atomName, aaName]=textread('InputFiles/TOCSY.m','%f %f %f %f %s %s');

fid = fopen('InputFiles/TOCSY.m.standardFormat','w');
for i = 1:length(TH1)
  fprintf(fid, '%f %f %d %d\n', TH1(i), TH2(i), TN(i), TRN(i));
end

fclose(fid);